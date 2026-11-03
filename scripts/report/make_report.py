#!/usr/bin/env python3
"""Build the MetaCheese run report (workflow step R.1, pipeline script 16).

Collects the results of every pipeline step from a run folder, joins the two
branches (community-level and assembly-specific profiling) and writes one
self-contained interactive HTML file (no internet connection or extra software
needed to open it), meant to be read in a browser and shared with colleagues.

Usage:
    python3 scripts/report/make_report.py output/<run> [--out FILE] [--samplesheet TSV]

Only the Python 3 standard library is used, so it runs both inside the
container and on the host. Steps whose outputs are missing are reported as
"not run" and the rest of the report is still produced.
"""
import argparse
import datetime as dt
import glob
import gzip
import html
import json
import math
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

# Tool versions installed in the metacheese image (see Dockerfile / README)
TOOLS = [
    ("Q.1", "FastQC", "0.12.1"),
    ("Q.1", "MultiQC", "1.25.1"),
    ("P.1", "AdapterRemoval", "2.3.3"),
    ("I.2, P.2, A.3", "Bowtie2", "2.3.3.1"),
    ("C.1", "MetaPhlAn", "4.1.0 (database mpa_vJun23_CHOCOPhlAnSGB_202307)"),
    ("C.2", "R + vegan / rbiom", "R 4.3.3"),
    ("C.3", "HUMAnN", "3.9 (ChocoPhlAn v201901_v31, UniRef90 DIAMOND)"),
    ("A.1", "metaSPAdes", "3.13.0"),
    ("A.3", "samtools", "1.6"),
    ("A.4", "MetaBAT2", "2.17"),
    ("A.5", "CheckM2", "1.0.2"),
    ("A.7", "TORMES", "1.3.0"),
    ("R.1", "make_report.py", "Python 3 standard library"),
]

# Pipeline scripts in execution order, with workflow IDs and output folders
STEPS = [
    ("00.fastqc", "Q.1", "Quality check of raw reads", ["00_fastqc"]),
    ("00-01.extract_and_adapter", "P.1", "Adapter and quality trimming", ["01_AdapterRemoval"]),
    ("03.bowtie-remove_host", "P.2", "Host DNA removal", ["03_bowtie2_output"]),
    ("04.metaphlan", "C.1, C.2", "Taxonomic profiling and diversity", ["04_metaphlan_output"]),
    ("04b-1.humann_job_array", "C.3", "Functional profiling", ["04b_humann_output"]),
    ("04b-last.humann_join_tables", "C.4", "Community-level aggregation", ["04b_humann_output/merged_tables"]),
    ("05.spades.assembly", "A.1", "Metagenome assembly", ["05_spades_output"]),
    ("06.contig_filter", "A.2", "Contig filtering", ["05_spades_output/contigs/filtered"]),
    ("07.bowtie-MAG_database", "A.3", "Contig index", ["07_Bowtie_Index"]),
    ("08-09.mapp_coverage", "A.3", "Read mapping and coverage", ["08_mapping_coverage"]),
    ("10-11-11b.metabat", "A.4", "Binning", ["11_metabat_MAG", "10_metabat_depth"]),
    ("12-13.checkm", "A.5", "MAG quality", ["13_checkm2"]),
    ("14.filter-metadata", "A.6", "MAG filtering and metadata", ["14_MAGs_high_quality"]),
    ("15.tormes_MAG", "A.7", "MAG annotation", ["15_tormes_MAGs"]),
    ("16.report", "R.1", "Integrated report", ["report"]),
]

RANKS = {"p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}


# --------------------------------------------------------------------------- helpers
def read_text(path, default=None):
    try:
        with open(path, errors="replace") as fh:
            return fh.read()
    except OSError:
        return default


def read_tsv(path, comment=None):
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or (comment and line.startswith(comment)):
                continue
            rows.append(line.split("\t"))
    return rows


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def clean_sample(name):
    """D_GP_02.unmapped.fq_Abundance-CPM -> D_GP_02 (sample names contain no dots;
    every pipeline file name adds its suffixes after a dot)"""
    return os.path.basename(name).split(".")[0]


def match_sample(col, samples):
    """Map a table column name to one of the known sample names."""
    base = clean_sample(col)
    if base in samples:
        return base
    for s in sorted(samples, key=len, reverse=True):
        if col.startswith(s + ".") or col.startswith(s + "_"):
            return s
    return base


def open_maybe_gz(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, errors="replace")


def fasta_lengths(path):
    lengths, n = [], 0
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lengths.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lengths.append(n)
    return lengths


def n50(lengths):
    if not lengths:
        return 0
    half, acc = sum(lengths) / 2, 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= half:
            return length
    return 0


def gc_content(path):
    gc = at = 0
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                s = line.strip().upper()
                gc += s.count("G") + s.count("C")
                at += s.count("A") + s.count("T")
    return 100.0 * gc / (gc + at) if gc + at else None


def jacobi_eigen(a, iters=100):
    """Eigen-decomposition of a small symmetric matrix (list of lists)."""
    n = len(a)
    a = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(iters):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-18:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-15:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    vals = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: -vals[i])
    return [vals[i] for i in order], [[v[k][i] for k in range(n)] for i in order]


def pcoa(dist):
    """Classical multidimensional scaling; returns 2D coordinates and % variance."""
    n = len(dist)
    if n < 3:
        return None
    d2 = [[dist[i][j] ** 2 for j in range(n)] for i in range(n)]
    rmean = [sum(r) / n for r in d2]
    gmean = sum(rmean) / n
    b = [[-0.5 * (d2[i][j] - rmean[i] - rmean[j] + gmean) for j in range(n)] for i in range(n)]
    vals, vecs = jacobi_eigen(b)
    pos = sum(x for x in vals if x > 0) or 1
    coords = [[vecs[k][i] * math.sqrt(max(vals[k], 0)) for k in range(2)] for i in range(n)]
    return {"coords": coords, "explained": [100 * max(vals[k], 0) / pos for k in range(2)]}


# --------------------------------------------------------------------------- collectors
def collect_samplesheet(run, config_text, explicit):
    candidates = [explicit] if explicit else []
    m = re.search(r'^step00:\s*\n\s*dir1:\s*"?([^"\n#]+)', config_text or "", re.M)
    if m:
        rel = m.group(1).strip()
        candidates.append(os.path.normpath(os.path.join(REPO, "scripts", "pipeline", rel, "samplesheet.tsv")))
        candidates.append(os.path.join(REPO, "input", os.path.basename(rel), "samplesheet.tsv"))
    for path in candidates:
        if path and os.path.isfile(path):
            rows = read_tsv(path)
            head = rows[0]
            out = []
            for r in rows[1:]:
                rec = dict(zip(head, r))
                out.append({"sample": rec.get("sample_name") or r[0],
                            "group": rec.get("group") or "",
                            "run": rec.get("run_accession", ""),
                            "title": rec.get("title", "")})
            return out, os.path.relpath(path, REPO) if path.startswith(REPO) else path
    return [], None


def collect_steps(run):
    times = {}
    path = os.path.join(run, "step_times.tsv")
    if os.path.isfile(path):
        for row in read_tsv(path):
            if len(row) >= 2 and to_float(row[1]) is not None:
                times[row[0].replace(".sh", "")] = float(row[1])
    out = []
    for script, sid, name, folders in STEPS:
        present = any(os.path.isdir(os.path.join(run, f)) and os.listdir(os.path.join(run, f)) for f in folders)
        out.append({"script": script, "id": sid, "name": name, "folder": folders[0],
                    "seconds": times.get(script), "present": present})
    return out


def collect_fastqc(run):
    res = []
    for z in sorted(glob.glob(os.path.join(run, "00_fastqc", "*_fastqc.zip"))):
        try:
            with zipfile.ZipFile(z) as zf:
                root = zf.namelist()[0].split("/")[0]
                summary = zf.read(f"{root}/summary.txt").decode()
                data = zf.read(f"{root}/fastqc_data.txt").decode()
        except (zipfile.BadZipFile, KeyError, IndexError):
            continue
        modules = [{"status": s, "module": m} for s, m, *_ in (l.split("\t") for l in summary.splitlines() if l)]
        basic = dict(re.findall(r"^(Total Sequences|Sequence length|%GC|Sequences flagged as poor quality)\t(.+)$", data, re.M))
        qual = []
        sec = re.search(r">>Per base sequence quality.*?\n#Base\tMean.*?\n(.*?)>>END_MODULE", data, re.S)
        if sec:
            for line in sec.group(1).strip().splitlines():
                f = line.split("\t")
                qual.append([f[0], round(float(f[1]), 2)])
        name = os.path.basename(z).replace("_fastqc.zip", "")
        mate = "R1" if re.search(r"_(1|L1)$", name) else "R2" if re.search(r"_(2|L2)$", name) else ""
        res.append({"file": name, "sample": re.sub(r"_(1|2|L1|L2)$", "", name), "mate": mate,
                    "total": int(basic.get("Total Sequences", 0) or 0),
                    "length": basic.get("Sequence length", ""), "gc": to_float(basic.get("%GC")),
                    "modules": modules, "quality": qual})
    return res


def collect_trimming(run):
    res = {}
    for f in sorted(glob.glob(os.path.join(run, "01_AdapterRemoval", "*.settings"))):
        t = read_text(f, "")
        get = lambda k: int(m.group(1)) if (m := re.search(rf"^{re.escape(k)}: (\d+)", t, re.M)) else None
        avg = re.search(r"^Average length of retained reads: ([\d.]+)", t, re.M)
        res[clean_sample(f.replace(".settings", ""))] = {
            "pairs": get("Total number of read pairs"),
            "adapters": sum(int(x) for x in re.findall(r"^Number of reads with adapters\[\d+\]: (\d+)", t, re.M)),
            "retained_reads": get("Number of retained reads"),
            "discarded": (get("Number of discarded mate 1 reads") or 0) + (get("Number of discarded mate 2 reads") or 0),
            "singletons": (get("Number of singleton mate 1 reads") or 0) + (get("Number of singleton mate 2 reads") or 0),
            "avg_length": float(avg.group(1)) if avg else None,
        }
    return res


def parse_bt2(text):
    m_in = re.search(r"^(\d+) reads; of these:", text, re.M)
    m_c0 = re.search(r"^\s*(\d+) \([\d.]+%\) aligned concordantly 0 times", text, re.M)
    m_rate = re.search(r"^([\d.]+)% overall alignment rate", text, re.M)
    if not m_in:
        return None
    return {"pairs_in": int(m_in.group(1)),
            "not_concordant": int(m_c0.group(1)) if m_c0 else None,
            "overall_rate": float(m_rate.group(1)) if m_rate else None}


def collect_host(run):
    res = {}
    for f in glob.glob(os.path.join(run, "03_bowtie2_output", "*.bowtie2.log")):
        r = parse_bt2(read_text(f, ""))
        if r:
            res[clean_sample(f)] = r
    if res:
        return res
    # Older runs: the Bowtie2 summaries are only in the run log
    for log in sorted(glob.glob(os.path.join(run, "run_*.log"))):
        text = read_text(log, "")
        for block in re.split(r"Starting mapping for ", text)[1:]:
            name = block.split()[0]
            r = parse_bt2(block.split("Executing")[0])
            if r and r["overall_rate"] is not None:
                res[name] = r
    return res


def collect_metaphlan(run, samples):
    out = {"profiles": {}, "reads": {}}
    path = os.path.join(run, "04_metaphlan_output", "merged_abundance_table.txt")
    if not os.path.isfile(path):
        return None
    rows = read_tsv(path, comment="#mpa")
    head = rows[0]
    cols = [match_sample(c, samples) for c in head[1:]]
    out["samples"] = cols
    for r in rows[1:]:
        clade = r[0].split("|")
        last = clade[-1]
        rank = RANKS.get(last[0])
        if not rank or last.startswith("t__"):
            continue
        name = last[3:].replace("_", " ")
        out["profiles"].setdefault(rank, []).append({"taxon": name, "lineage": " › ".join(c[3:].replace("_", " ") for c in clade[:-1]),
                                                     "values": [to_float(x) or 0.0 for x in r[1:]]})
    for f in glob.glob(os.path.join(run, "04_metaphlan_output", "*.txt")):
        m = re.search(r"^#(\d+) reads processed", read_text(f, ""), re.M)
        if m:
            out["reads"][clean_sample(f)] = int(m.group(1))
    return out


def collect_diversity(run, samples):
    d = os.path.join(run, "04_metaphlan_output", "diversity")
    if not os.path.isdir(d):
        return None
    alpha, beta = {}, {}
    labels = {"richness": "Observed species", "shannon": "Shannon index", "simpson": "Gini-Simpson index", "gini": "Gini dominance"}
    for key, label in labels.items():
        f = os.path.join(d, f"merged_abundance_table_{key}.tsv")
        if os.path.isfile(f):
            rows = read_tsv(f)
            alpha[key] = {"label": label, "values": {match_sample(r[0], samples): to_float(r[1]) for r in rows[1:] if len(r) > 1}}
    for key, label in (("bray-curtis", "Bray-Curtis"), ("jaccard", "Jaccard")):
        f = os.path.join(d, f"merged_abundance_table_{key}.tsv")
        if os.path.isfile(f):
            rows = read_tsv(f)
            names = [match_sample(c, samples) for c in rows[0]]
            mat = [[to_float(x) or 0.0 for x in r[1:]] for r in rows[1:]]
            beta[key] = {"label": label, "samples": names, "matrix": mat, "pcoa": pcoa(mat)}
    return {"alpha": alpha, "beta": beta}


def read_humann_table(path, samples, top=None):
    rows = read_tsv(path)
    if not rows:
        return None
    cols = [match_sample(c, samples) for c in rows[0][1:]]
    feats = []
    special = {}
    for r in rows[1:]:
        vals = [to_float(x) or 0.0 for x in r[1:]]
        if r[0] in ("UNMAPPED", "UNINTEGRATED", "UNGROUPED"):
            special[r[0]] = vals
            continue
        if "|" in r[0]:
            continue
        name = r[0]
        ident, _, desc = name.partition(": ")
        feats.append({"id": ident, "name": desc or ident, "values": vals})
    total = len(feats)
    feats.sort(key=lambda f: -sum(f["values"]))
    return {"samples": cols, "features": feats[:top] if top else feats, "total": total, "special": special}


def humann_sample_stats(path, special):
    """Unstratified features of one per-sample HUMAnN table: count of real
    features and share (%) of each special row (UNMAPPED, UNINTEGRATED)."""
    total, n, sp = 0.0, 0, {k: 0.0 for k in special}
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or "|" in line.split("\t", 1)[0]:
                continue
            f = line.rstrip("\n").split("\t")
            v = to_float(f[1]) if len(f) > 1 else None
            if v is None:
                continue
            total += v
            if f[0] in sp:
                sp[f[0]] += v
            else:
                n += 1
    return n, {k: (100 * v / total if total else None) for k, v in sp.items()}


def humann_species_contribution(path, samples, top=10):
    """Share (%) of the stratified pathway abundance carried by each species,
    per sample: which organisms the community's functions come from."""
    rows = read_tsv(path)
    if not rows:
        return None
    cols = [match_sample(c, samples) for c in rows[0][1:]]
    per = {}
    for r in rows[1:]:
        feat, _, strat = r[0].partition("|")
        if not strat or feat in ("UNMAPPED", "UNINTEGRATED"):
            continue
        sp = strat.split(".s__")[-1] if ".s__" in strat else strat
        sp = sp.replace("_", " ")
        acc = per.setdefault(sp, [0.0] * len(cols))
        for i, x in enumerate(r[1:]):
            acc[i] += to_float(x) or 0.0
    tot = [sum(v[i] for v in per.values()) for i in range(len(cols))]
    share = {k: [100 * v[i] / tot[i] if tot[i] else 0.0 for i in range(len(cols))] for k, v in per.items()}
    order = sorted(share, key=lambda k: -sum(share[k]))
    return {"samples": cols, "species": [{"taxon": k, "values": share[k]} for k in order[:top]], "total": len(order)}


def collect_humann(run, samples):
    base = os.path.join(run, "04b_humann_output")
    if not os.path.isdir(base):
        return None
    out = {"per_sample": sorted(clean_sample(f) for f in glob.glob(os.path.join(base, "*_genefamilies.tsv"))), "tables": {},
           "stats": {}, "contrib": None}
    # C.3: per-sample results of step 04b-1
    for f in glob.glob(os.path.join(base, "*_genefamilies.tsv")):
        s = match_sample(os.path.basename(f), samples)
        n_gf, sp_gf = humann_sample_stats(f, ["UNMAPPED"])
        st = {"genefamilies": n_gf, "unmapped": sp_gf["UNMAPPED"]}
        pw = f.replace("_genefamilies.tsv", "_pathabundance.tsv")
        if os.path.isfile(pw):
            n_pw, sp_pw = humann_sample_stats(pw, ["UNMAPPED", "UNINTEGRATED"])
            st.update({"pathways": n_pw, "unintegrated": sp_pw["UNINTEGRATED"]})
        out["stats"][s] = st
    # C.4: merged tables of step 04b-last
    strat = os.path.join(base, "merged_tables", "humann_joined_pathabundance_cpm.tsv")
    if os.path.isfile(strat):
        out["contrib"] = humann_species_contribution(strat, samples)
    un = os.path.join(base, "merged_tables", "unstratified")
    specs = [("pathways", "MetaCyc pathways", "humann_joined_pathabundance_cpm.tsv", "CPM"),
             ("ko", "KEGG Orthology groups", "merged_table_ko.tsv", "CPM"),
             ("ec", "Enzymes (EC numbers)", "merged_table_ec.tsv", "CPM"),
             ("genefamilies", "UniRef90 gene families", "humann_joined_genefamilies_cpm.tsv", "CPM")]
    for key, label, fname, unit in specs:
        f = os.path.join(un, fname)
        if os.path.isfile(f):
            t = read_humann_table(f, samples, top=0 if key == "genefamilies" else 40)
            if key == "genefamilies":
                t["features"] = []
            t.update({"label": label, "unit": unit})
            out["tables"][key] = t
    return out


def collect_assembly(run):
    res = {}
    for f in sorted(glob.glob(os.path.join(run, "05_spades_output", "contigs", "filtered", "*.fasta"))):
        lengths = fasta_lengths(f)
        raw = os.path.join(run, "05_spades_output", "contigs", os.path.basename(f).replace("_sort.fasta", ""))
        raw_lengths = fasta_lengths(raw) if os.path.isfile(raw) else []
        res[clean_sample(f)] = {"contigs": len(lengths), "total_bp": sum(lengths), "n50": n50(lengths),
                                "longest": max(lengths) if lengths else 0, "raw_contigs": len(raw_lengths),
                                "raw_bp": sum(raw_lengths)}
    return res


def collect_mapping(run, samples):
    res = {}
    # step 08-09 writes <sample>_bowtie2.log (no dot), so match on the prefix
    for f in glob.glob(os.path.join(run, "08_mapping_coverage", "*_bowtie2.log")):
        r = parse_bt2(read_text(f, ""))
        if r:
            res[match_sample(os.path.basename(f), samples)] = r["overall_rate"]
    return res


def collect_mags(run, samples):
    qual = {}
    f = os.path.join(run, "13_checkm2", "quality_report.tsv")
    if os.path.isfile(f):
        rows = read_tsv(f)
        head = rows[0]
        for r in rows[1:]:
            rec = dict(zip(head, r))
            qual[rec["Name"]] = rec
    hq_dir = os.path.join(run, "14_MAGs_high_quality")
    hq = {os.path.basename(p)[:-3] for p in glob.glob(os.path.join(hq_dir, "*.fa"))}
    files = {os.path.basename(p)[:-3]: p for p in glob.glob(os.path.join(run, "11_metabat_MAG", "*.fa"))}
    files.update({os.path.basename(p)[:-3]: p for p in glob.glob(os.path.join(hq_dir, "*.fa"))})
    names = sorted(set(qual) | set(files))
    mags = []
    for n in names:
        q = qual.get(n, {})
        size = to_float(q.get("Genome_Size"))
        if size is None and n in files:
            size = sum(fasta_lengths(files[n]))
        mags.append({"name": n, "sample": match_sample(n, samples),
                     "completeness": to_float(q.get("Completeness")),
                     "contamination": to_float(q.get("Contamination")),
                     "size": size, "n50": to_float(q.get("Contig_N50")),
                     "gc": to_float(q.get("GC_Content")) if q.get("GC_Content") else None,
                     "model": q.get("Completeness_Model_Used", ""), "passed": n in hq})
    return mags


def collect_tormes(run):
    d = os.path.join(run, "15_tormes_MAGs")
    if not os.path.isdir(d):
        return None
    reports = sorted(glob.glob(os.path.join(d, "**", "*.html"), recursive=True))
    return {"reports": [os.path.relpath(p, run) for p in reports][:20]}


def fasta_ids(path):
    with open_maybe_gz(path) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def collect_tormes_taxonomy(run):
    """MAG -> (family, genus, species) from TORMES' Kraken2 identification."""
    tax = {}
    for f in glob.glob(os.path.join(run, "15_tormes_MAGs", "**", "taxonomic-identification-kraken2.txt"), recursive=True):
        rows = read_tsv(f)
        for r in rows[1:]:
            if len(r) < 4:
                continue
            strip = lambda x: re.sub(r"\s*\(.*?\)\s*$", "", x).strip()
            tax[r[0]] = {"family": strip(r[1]), "genus": strip(r[2]), "species": strip(r[3])}
    return tax


def collect_integration(run, samples, mags, metaphlan, mapping):
    """R.1: join the two branches.

    For each sample, the reads that mapped to the sample's own contigs (A.3)
    are split by where their contigs ended up: kept MAGs (A.6), other bins
    (A.4), contigs that were not binned. Contig depth x length (A.3 depth
    file) approximates the mapped bases. Each MAG gets its share of the
    microbial reads and, when TORMES classified it (A.7), the abundance that
    MetaPhlAn (C.1) gives to the same species or genus in the same sample.
    """
    bins = {}
    for d, kept in (("11_metabat_MAG", False), ("14_MAGs_high_quality", True)):
        for f in glob.glob(os.path.join(run, d, "*.fa")):
            name = os.path.basename(f)[:-3]
            s = match_sample(name, samples)
            for c in fasta_ids(f):
                bins[(s, c)] = (name, kept)  # SPAdes contig names repeat across samples
    fate, mag_share = {}, {}
    for f in glob.glob(os.path.join(run, "10_metabat_depth", "*.depth.txt")):
        s = match_sample(os.path.basename(f), samples)
        rows = read_tsv(f)
        if not rows:
            continue
        tot = kept = other = 0.0
        per_bin = {}
        for r in rows[1:]:
            ln, dp = to_float(r[1]), to_float(r[2])
            if ln is None or dp is None:
                continue
            b = ln * dp
            tot += b
            hit = bins.get((s, r[0]))
            if hit:
                per_bin[hit[0]] = per_bin.get(hit[0], 0.0) + b
                if hit[1]:
                    kept += b
                else:
                    other += b
        if not tot:
            continue
        mapped = mapping.get(s)
        scale = (mapped if mapped is not None else 100.0) / tot
        fate[s] = {"kept": kept * scale, "other_bins": other * scale, "unbinned": (tot - kept - other) * scale,
                   "not_assembled": 100.0 - mapped if mapped is not None else None}
        for n, b in per_bin.items():
            mag_share[n] = b * scale

    tax = collect_tormes_taxonomy(run)
    prof = {}
    if metaphlan:
        cols = metaphlan["samples"]
        for rank in ("species", "genus"):
            for t in metaphlan["profiles"].get(rank, []):
                for i, x in enumerate(cols):
                    prof[(rank, x, t["taxon"].lower())] = t["values"][i]
    norm = lambda x: re.sub(r"\s+", " ", (x or "").replace("_", " ")).strip().lower()
    out_mags = []
    for m in mags:
        t = tax.get(m["name"], {})
        sp, gn = norm(t.get("species")), norm(t.get("genus"))
        mpa_sp = prof.get(("species", m["sample"], sp)) if sp else None
        mpa_gn = prof.get(("genus", m["sample"], gn)) if gn else None
        if mpa_gn is None and sp:
            mpa_gn = prof.get(("genus", m["sample"], sp.split(" ")[0]))
        out_mags.append({"name": m["name"], "sample": m["sample"], "passed": m["passed"],
                         "completeness": m["completeness"], "contamination": m["contamination"],
                         "share": mag_share.get(m["name"]), "species": t.get("species"), "genus": t.get("genus"),
                         "mpa_species": mpa_sp, "mpa_genus": mpa_gn})
    return {"fate": fate, "mags": out_mags, "taxonomy": bool(tax)}


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="run folder, e.g. output/example_small")
    ap.add_argument("--out", help="output HTML (default <run>/report/report.html)")
    ap.add_argument("--samplesheet", help="TSV with sample_name and group columns (default: samplesheet.tsv in the input folder)")
    ap.add_argument("--title", help="report title (default: run name)")
    ap.add_argument("--fragment", action="store_true", help="omit <html>/<head>/<body> (for embedding in another page)")
    args = ap.parse_args()

    run = os.path.abspath(args.run)
    if not os.path.isdir(run):
        sys.exit(f"ERROR: run folder {run} not found")
    name = os.path.basename(run)
    config_text = read_text(os.path.join(run, "config.yml"), "")

    sheet, sheet_path = collect_samplesheet(run, config_text, args.samplesheet)
    trimming = collect_trimming(run)
    host = collect_host(run)
    samples = [s["sample"] for s in sheet] or sorted(set(trimming) | set(host))
    groups = {s["sample"]: s["group"] for s in sheet}

    metaphlan = collect_metaphlan(run, samples)
    mapping = collect_mapping(run, samples)
    mags = collect_mags(run, samples)
    data = {
        "run": name,
        "title": args.title or name,
        "generated": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "samplesheet": sheet_path,
        "samples": [{"sample": s, "group": groups.get(s, ""),
                     "run": next((x["run"] for x in sheet if x["sample"] == s), ""),
                     "title": next((x["title"] for x in sheet if x["sample"] == s), "")} for s in samples],
        "steps": collect_steps(run),
        "fastqc": collect_fastqc(run),
        "trimming": trimming,
        "host": host,
        "metaphlan": metaphlan,
        "diversity": collect_diversity(run, samples),
        "humann": collect_humann(run, samples),
        "assembly": collect_assembly(run),
        "mapping": mapping,
        "mags": mags,
        "integration": collect_integration(run, samples, mags, metaphlan, mapping),
        "tormes": collect_tormes(run),
        "tools": [{"id": a, "tool": b, "version": c} for a, b, c in TOOLS],
        "config": config_text,
    }

    svg = read_text(os.path.join(REPO, "docs", "img", "workflow.svg"), "")
    svg = re.sub(r'\s(width|height)="\d+"', "", svg, count=2)
    template = read_text(os.path.join(HERE, "report_template.html"))
    if template is None:
        sys.exit("ERROR: report_template.html not found next to make_report.py")
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    page = (template.replace("/*__DATA__*/null", payload)
                    .replace("<!--__WORKFLOW__-->", svg)
                    .replace("__TITLE__", html.escape(data["title"])))
    if args.fragment:
        page = "".join(re.findall(r"(?s)<!--__FRAGMENT_START__-->(.*?)<!--__FRAGMENT_END__-->", page))
    else:
        page = page.replace("<!--__FRAGMENT_START__-->", "").replace("<!--__FRAGMENT_END__-->", "")
    out = args.out or os.path.join(run, "report", "report.html")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as fh:
        fh.write(page)
    print(f"Report written to {out} ({os.path.getsize(out) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()

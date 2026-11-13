#!/usr/bin/env python3
"""Joint analyses of the two branches (workflow step R.1, pipeline script 16).

Combines the results of the community-level branch (C.1-C.4) and of the
assembly-specific branch (A.1-A.7) of one run and writes them as plain tables
in <run>/16_integration/, read by the final report (R.2) and usable on their own:

    read_fate.tsv  one row per sample: where the sample's microbial reads end up
                   (kept MAGs, other bins, unbinned contigs, not assembled), as %
    mags.tsv       one row per MAG (bin): quality, share of the sample's
                   microbial reads, Kraken2 taxonomy (TORMES), and the abundance
                   MetaPhlAn gives to the same species and genus in the same sample
    notes.txt      notes on how names were matched (one per line, may be empty)

Functions of the kept MAGs (A.7, Prokka annotation by TORMES) against the
functional profiles (C.3, HUMAnN), linked in two ways:

    by species     HUMAnN splits each sample's pathways and gene families by the
                   species they come from; a MAG is linked to the HUMAnN species
                   with the same name as its Kraken2 species (both databases use
                   the names from before the 2020 split of Lactobacillus)
    by EC number   the EC numbers Prokka finds in the MAG against those HUMAnN
                   finds in the sample, for the MAG's species and for the whole
                   community (C.4 table merged_table_ec.tsv, which keeps the
                   species split of C.3)

    mag_functions.tsv  one row per kept MAG: the linked HUMAnN species, its
                       pathways and gene families and their share of the sample,
                       and EC counts (in the MAG, in HUMAnN's species, in both,
                       in only one, in the community)
    mag_pathways.tsv   one row per MAG and pathway HUMAnN gives to its species:
                       abundance, and share of the pathway's community abundance
    mag_ecs.tsv        one row per MAG and EC number found in the MAG or in
                       HUMAnN's species: in the MAG or not, CPM in the species
                       and in the community
    ec_coverage.tsv    one row per sample: EC numbers of the community, how many
                       the kept MAGs encode, and the share of the community's EC
                       abundance they encode

Usage:
    python3 scripts/integration/integrate.py output/<run> [--samplesheet TSV]

Only the Python 3 standard library is used. The parsers of the other steps'
outputs are those of the report (scripts/report/make_report.py), so both read
the files the same way. Inputs that are missing give empty tables.
"""
import argparse
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "report"))
import make_report as mr  # noqa: E402

OUT_DIR = "16_integration"

# Genera that reference databases of different ages name differently. A MAG whose
# TORMES genus is a key here also matches a detected species of the same epithet
# in any genus accepted by "current". If no species matches, its genus abundance
# is the sum of those genera, because the old name cannot say which new genus the
# MAG belongs to. Add an entry when another split matters.
RENAMED_GENERA = {
    "lactobacillus": {
        "current": re.compile(r"^(\w*lactobacillus|lacti\w*bacillus)$"),
        "note": "TORMES' Kraken2 database uses the names from before the 2020 split of Lactobacillus "
                "(Zheng et al. 2020; e.g. Lactobacillus rhamnosus), while MetaPhlAn uses the new genera "
                "(Lacticaseibacillus rhamnosus). These MAGs are matched by species epithet across the "
                "split, and the MetaPhlAn name is shown when it differs.",
        "genus_note": "For a Lactobacillus MAG with no MetaPhlAn name, the MetaPhlAn genus value is the sum of "
                      "Lactobacillus and the genera split from it, because the old name cannot say which new "
                      "genus the MAG belongs to.",
    },
}

FATE_COLUMNS = [("kept", "kept_mags_pct"), ("other_bins", "other_bins_pct"),
                ("unbinned", "unbinned_pct"), ("not_assembled", "not_assembled_pct")]
MAG_COLUMNS = [("name", "mag"), ("sample", "sample"), ("passed", "kept"),
               ("completeness", "completeness"), ("contamination", "contamination"),
               ("share", "read_share_pct"), ("species", "kraken2_species"), ("genus", "kraken2_genus"),
               ("mpa_name", "metaphlan_name"), ("mpa_species", "metaphlan_species_pct"),
               ("mpa_genus", "metaphlan_genus_pct")]


def read_fate(run, samples, mapping):
    """For each sample, the reads that mapped to the sample's own contigs (A.3)
    are split by where their contigs ended up: kept MAGs (A.6), other bins
    (A.4), contigs that were not binned. Contig depth x length (A.4 depth
    file) approximates the mapped bases. Also returns each bin's share of the
    sample's microbial reads."""
    bins = {}
    for d, kept in (("11_metabat_MAG", False), ("14_MAGs_high_quality", True)):
        for f in glob.glob(os.path.join(run, d, "*.fa")):
            name = os.path.basename(f)[:-3]
            s = mr.match_sample(name, samples)
            for c in mr.fasta_ids(f):
                bins[(s, c)] = (name, kept)  # SPAdes contig names repeat across samples
    fate, mag_share = {}, {}
    for f in glob.glob(os.path.join(run, "10_metabat_depth", "*.depth.txt")):
        s = mr.match_sample(os.path.basename(f), samples)
        rows = mr.read_tsv(f)
        if not rows:
            continue
        tot = kept = other = 0.0
        per_bin = {}
        for r in rows[1:]:
            ln, dp = mr.to_float(r[1]), mr.to_float(r[2])
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
    return fate, mag_share


def mags_vs_metaphlan(run, mags, metaphlan, mag_share):
    """Each MAG with its read share and, when TORMES classified it (A.7), the
    abundance that MetaPhlAn (C.1) gives to the same species or genus in the
    same sample."""
    tax = mr.collect_tormes_taxonomy(run)
    prof = {}
    if metaphlan:
        cols = metaphlan["samples"]
        for rank in ("species", "genus"):
            for t in metaphlan["profiles"].get(rank, []):
                for i, x in enumerate(cols):
                    # The merged table lists every taxon in every sample: keep
                    # only those MetaPhlAn really finds in this sample.
                    if (t["values"][i] or 0) > 0:
                        prof[(rank, x, t["taxon"].lower())] = t["values"][i]
    norm = lambda x: re.sub(r"\s+", " ", (x or "").replace("_", " ")).strip().lower()
    by_epithet = {}
    for (rank, x, taxon) in prof:
        w = taxon.split(" ")
        if rank == "species" and len(w) == 2:
            for old, rule in RENAMED_GENERA.items():
                if rule["current"].match(w[0]):
                    by_epithet[(old, x, w[1])] = taxon
    notes = []
    out = []
    for m in mags:
        t = tax.get(m["name"], {})
        sp, gn = norm(t.get("species")), norm(t.get("genus"))
        mpa_name = sp if sp and ("species", m["sample"], sp) in prof else None
        if not mpa_name and gn in RENAMED_GENERA and len(sp.split(" ")) == 2:
            mpa_name = by_epithet.get((gn, m["sample"], sp.split(" ")[1]))
            note = RENAMED_GENERA[gn]["note"]
            if mpa_name and note not in notes:
                notes.append(note)
        mpa_sp = prof.get(("species", m["sample"], mpa_name)) if mpa_name else None
        mpa_gn = None
        if mpa_name:
            mpa_gn = prof.get(("genus", m["sample"], mpa_name.split(" ")[0]))
        if mpa_gn is None and gn in RENAMED_GENERA:
            rule = RENAMED_GENERA[gn]
            parts = [v for (rank, x, taxon), v in prof.items()
                     if rank == "genus" and x == m["sample"] and rule["current"].match(taxon)]
            if parts:
                mpa_gn = sum(parts)
                if rule["genus_note"] not in notes:
                    notes.append(rule["genus_note"])
        elif mpa_gn is None and gn:
            mpa_gn = prof.get(("genus", m["sample"], gn))
        out.append({"name": m["name"], "sample": m["sample"], "passed": m["passed"],
                    "completeness": m["completeness"], "contamination": m["contamination"],
                    "share": mag_share.get(m["name"]), "species": t.get("species"), "genus": t.get("genus"),
                    "mpa_name": mpa_name[0].upper() + mpa_name[1:] if mpa_name else None,
                    "mpa_species": mpa_sp, "mpa_genus": mpa_gn})
    return out, notes


SPECIAL_ROWS = ("UNMAPPED", "UNINTEGRATED", "UNGROUPED")
FULL_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")  # partial ECs (1.1.1.-) cannot be matched


def species_key(name):
    """Comparable species name: 'g__Lactobacillus.s__Lactobacillus_casei' and
    'Lactobacillus casei' both give 'lactobacillus casei'."""
    name = (name or "").split(".s__")[-1]
    return re.sub(r"\s+", " ", name.replace("_", " ")).strip().lower()


def species_label(strat):
    """HUMAnN stratum to a display name; None for 'unclassified'."""
    if not strat or strat == "unclassified":
        return None
    return strat.split(".s__")[-1].replace("_", " ")


def humann_by_species(run, samples):
    """C.3 per-sample tables split by species. For each sample: pathway
    abundances per species and for the community, and gene-family counts and
    abundances per species."""
    out = {}
    base = os.path.join(run, "04b_humann_output")
    for f in glob.glob(os.path.join(base, "*_pathabundance.tsv")):
        s = mr.match_sample(os.path.basename(f), samples)
        rec = out.setdefault(s, {})
        community, per, names, total = {}, {}, {}, 0.0
        for r in mr.read_tsv(f, comment="#"):
            if len(r) < 2:
                continue
            feat, _, strat = r[0].partition("|")
            if feat in SPECIAL_ROWS:
                continue
            v = mr.to_float(r[1]) or 0.0
            ident, _, desc = feat.partition(": ")
            names[ident] = desc or ident
            if not strat:
                community[ident] = v
                continue
            total += v
            label = species_label(strat)
            if label and v > 0:
                per.setdefault(species_key(label), {"label": label, "values": {}})["values"][ident] = v
        rec.update({"pathways": per, "pathway_community": community, "pathway_names": names,
                    "pathway_stratified_total": total})
    for f in glob.glob(os.path.join(base, "*_genefamilies.tsv")):
        s = mr.match_sample(os.path.basename(f), samples)
        rec = out.setdefault(s, {})
        per, total = {}, 0.0
        with open(f, errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                r = line.rstrip("\n").split("\t")
                if len(r) < 2:
                    continue
                v = mr.to_float(r[1]) or 0.0
                feat, _, strat = r[0].partition("|")
                if not strat:
                    total += v  # all reads, UNMAPPED included
                    continue
                label = species_label(strat)
                if label and v > 0:
                    acc = per.setdefault(species_key(label), [0, 0.0])
                    acc[0] += 1
                    acc[1] += v
        rec.update({"genefamilies": per, "genefamily_total": total})
    return out


def humann_ecs(run, samples):
    """EC numbers per sample (CPM) from the C.4 regrouped table, which keeps
    HUMAnN's species split: community values, values per species, names."""
    f = os.path.join(run, "04b_humann_output", "merged_tables", "merged_table_ec.tsv")
    if not os.path.isfile(f):
        return None
    rows = mr.read_tsv(f)
    cols = [mr.match_sample(c, samples) for c in rows[0][1:]]
    community = {s: {} for s in cols}
    per = {s: {} for s in cols}
    names = {}
    for r in rows[1:]:
        feat, _, strat = r[0].partition("|")
        if feat in SPECIAL_ROWS:
            continue
        ec, _, desc = feat.partition(": ")
        names[ec] = desc or ec
        label = species_label(strat) if strat else None
        for s, x in zip(cols, r[1:]):
            v = mr.to_float(x) or 0.0
            if v <= 0:
                continue
            if not strat:
                community[s][ec] = v
            elif label:
                per[s].setdefault(species_key(label), {})[ec] = v
    return {"community": community, "species": per, "names": names}


def mag_ecs(run):
    """Full EC numbers of each MAG from TORMES' Prokka tables."""
    out = {}
    for f in glob.glob(os.path.join(run, "15_tormes_MAGs", "annotation", "*_annotation", "*.tsv")):
        rows = mr.read_tsv(f)
        if not rows or "EC_number" not in rows[0]:
            continue
        i = rows[0].index("EC_number")
        ecs = set()
        for r in rows[1:]:
            if len(r) > i:
                ecs.update(e for e in re.split(r"[\s,;]+", r[i]) if FULL_EC.match(e))
        out[os.path.basename(f)[:-4]] = ecs
    return out


def mag_functions(run, samples, joint):
    """Kept MAGs against the functional profiles: by species and by EC number.
    Returns the rows of mag_functions, mag_pathways, mag_ecs and ec_coverage."""
    hum = humann_by_species(run, samples)
    ec = humann_ecs(run, samples)
    annotated = mag_ecs(run)
    summary, pathways, ecs = [], [], []
    kept = [m for m in joint if m["passed"]]
    for m in kept:
        s, key = m["sample"], species_key(m["species"])
        h = hum.get(s, {})
        pw = h.get("pathways", {}).get(key) if key else None
        gf = h.get("genefamilies", {}).get(key) if key else None
        sp_ecs = (ec["species"].get(s, {}).get(key) if ec and key else None)
        linked = pw is not None or gf is not None or sp_ecs is not None
        own = annotated.get(m["name"])
        community = ec["community"].get(s, {}) if ec else {}
        row = {"mag": m["name"], "sample": s, "kraken2_species": m["species"],
               "humann_species": (pw or {}).get("label") or (species_label(m["species"]) if linked else None)}
        if pw is not None:
            tot = h.get("pathway_stratified_total") or 0.0
            row["humann_pathways"] = len(pw["values"])
            row["humann_pathway_share_pct"] = 100 * sum(pw["values"].values()) / tot if tot else None
            com = h.get("pathway_community", {})
            for ident, v in sorted(pw["values"].items(), key=lambda kv: -kv[1]):
                pathways.append([m["name"], s, ident, h["pathway_names"].get(ident, ident), v,
                                 100 * v / com[ident] if com.get(ident) else None])
        elif linked:
            row["humann_pathways"] = 0
        if gf is not None:
            tot = h.get("genefamily_total") or 0.0
            row["humann_gene_families"] = gf[0]
            row["humann_gene_family_share_pct"] = 100 * gf[1] / tot if tot else None
        elif linked:
            row["humann_gene_families"] = 0
        if own is not None:
            row["mag_ecs"] = len(own)
            if ec:
                row["mag_ecs_in_community"] = len(own & set(community))
        if sp_ecs is not None or (linked and ec):
            sp = set(sp_ecs or {})
            row["humann_species_ecs"] = len(sp)
            if own is not None:
                row.update({"ecs_in_both": len(own & sp), "ecs_humann_species_only": len(sp - own),
                            "ecs_mag_only": len(own - sp)})
        if own is not None or sp_ecs:
            for e in sorted((own or set()) | set(sp_ecs or {})):
                ecs.append([m["name"], s, e, (ec["names"].get(e) if ec else None) or "",
                            None if own is None else (e in own),
                            (sp_ecs or {}).get(e, 0.0) if (linked and ec) else None,
                            community.get(e, 0.0) if ec else None])
        summary.append(row)
    coverage = []
    if ec and annotated:
        for s in [x for x in samples if x in ec["community"]] + sorted(set(ec["community"]) - set(samples)):
            com = ec["community"][s]
            encoded = set().union(*[annotated.get(m["name"], set()) for m in kept if m["sample"] == s])
            tot = sum(com.values())
            coverage.append([s, len(com), len(set(com) & encoded),
                             100 * sum(v for e, v in com.items() if e in encoded) / tot if tot else None])
    return summary, pathways, ecs, coverage


FUNCTION_COLUMNS = ["mag", "sample", "kraken2_species", "humann_species", "humann_pathways",
                    "humann_pathway_share_pct", "humann_gene_families", "humann_gene_family_share_pct",
                    "mag_ecs", "humann_species_ecs", "ecs_in_both", "ecs_humann_species_only", "ecs_mag_only",
                    "mag_ecs_in_community"]


def cell(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_tsv(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(cell(v) for v in r) + "\n")
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="run folder, e.g. output/example_small")
    ap.add_argument("--samplesheet", help="TSV with a sample_name column (default: samplesheet.tsv in the input folder)")
    args = ap.parse_args()

    run = os.path.abspath(args.run)
    if not os.path.isdir(run):
        sys.exit(f"ERROR: run folder {run} not found")
    config_text = mr.read_text(os.path.join(run, "config.yml"), "")
    sheet, _ = mr.collect_samplesheet(run, config_text, args.samplesheet)
    samples = [s["sample"] for s in sheet] or sorted(set(mr.collect_trimming(run)) | set(mr.collect_host(run)))

    metaphlan = mr.collect_metaphlan(run, samples)
    mapping = mr.collect_mapping(run, samples)
    mags = mr.collect_mags(run, samples)
    fate, mag_share = read_fate(run, samples, mapping)
    joint, notes = mags_vs_metaphlan(run, mags, metaphlan, mag_share)

    out = os.path.join(run, OUT_DIR)
    os.makedirs(out, exist_ok=True)
    order = [s for s in samples if s in fate] + sorted(s for s in fate if s not in samples)
    write_tsv(os.path.join(out, "read_fate.tsv"), ["sample"] + [c for _, c in FATE_COLUMNS],
              [[s] + [fate[s][k] for k, _ in FATE_COLUMNS] for s in order])
    write_tsv(os.path.join(out, "mags.tsv"), [c for _, c in MAG_COLUMNS],
              [[m[k] for k, _ in MAG_COLUMNS] for m in joint])
    with open(os.path.join(out, "notes.txt"), "w") as fh:
        fh.write("".join(n + "\n" for n in notes))
    functions, pathways, ecs, coverage = mag_functions(run, samples, joint)
    write_tsv(os.path.join(out, "mag_functions.tsv"), FUNCTION_COLUMNS,
              [[r.get(c) for c in FUNCTION_COLUMNS] for r in functions])
    write_tsv(os.path.join(out, "mag_pathways.tsv"),
              ["mag", "sample", "pathway", "pathway_name", "abundance", "share_of_community_pct"], pathways)
    write_tsv(os.path.join(out, "mag_ecs.tsv"),
              ["mag", "sample", "ec", "ec_name", "in_mag", "humann_species_cpm", "humann_community_cpm"], ecs)
    write_tsv(os.path.join(out, "ec_coverage.tsv"),
              ["sample", "community_ecs", "community_ecs_in_kept_mags", "ec_abundance_in_kept_mags_pct"], coverage)
    print(f"Joint analyses written to {out}: {len(fate)} samples, {len(joint)} MAGs "
          f"({sum(m['passed'] for m in joint)} kept, {sum(m['mpa_species'] is not None for m in joint if m['passed'])} "
          f"of them matching a MetaPhlAn species); {sum(r.get('humann_species') is not None for r in functions)} "
          f"kept MAGs linked to a HUMAnN species")


if __name__ == "__main__":
    main()

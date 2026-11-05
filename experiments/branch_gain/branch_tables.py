#!/usr/bin/env python3
"""Measure what each MetaCheese branch would miss on its own, for one run.

Uses only files in the run folder, and the same joins as the run report
(scripts/report/make_report.py). Writes these tables to <outdir>:

  mag_vs_metaphlan.tsv     kept MAGs: read share, share among kept MAGs, MetaPhlAn
                           abundance of the species TORMES names (0 = not found)
  species_without_mag.tsv  MetaPhlAn species >= --min-abundance with no kept MAG
  dominant_species.tsv     top species per sample according to each branch
  read_fate.tsv            microbial reads per sample: kept MAGs, other bins,
                           unbinned, not assembled; share on high-copy contigs
  humann_unmapped.tsv      HUMAnN reads with no gene family, per sample
  summary.txt              the headline numbers

Usage: python3 branch_tables.py <run_folder> <outdir> [--min-abundance 1] [--high-copy-depth 2000]
"""
import argparse
import glob
import importlib.util
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "make_report", os.path.join(HERE, "..", "..", "scripts", "report", "make_report.py"))
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)


def write(path, head, rows):
    with open(path, "w") as o:
        o.write("\t".join(head) + "\n")
        for r in rows:
            o.write("\t".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + "\n")


def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / den if den else float("nan")


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            out[order[k]] = (i + j) / 2
        i = j + 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("outdir")
    ap.add_argument("--min-abundance", type=float, default=1.0)
    ap.add_argument("--high-copy-depth", type=float, default=2000.0)
    a = ap.parse_args()
    run = os.path.abspath(a.run)
    os.makedirs(a.outdir, exist_ok=True)

    # Samples come from the run itself: the input folder may have changed since.
    samples = sorted(os.path.basename(f).split(".")[0]
                     for f in glob.glob(os.path.join(run, "10_metabat_depth", "*.depth.txt")))
    if not samples:
        sys.exit(f"ERROR: no depth files in {run}/10_metabat_depth")
    metaphlan = mr.collect_metaphlan(run, samples)
    mapping = mr.collect_mapping(run, samples)
    mags = mr.collect_mags(run, samples)
    ig = mr.collect_integration(run, samples, mags, metaphlan, mapping)

    kept = [m for m in ig["mags"] if m["passed"]]
    rows, x, y = [], [], []
    for s in samples:
        ms = [m for m in kept if m["sample"] == s]
        total = sum(m["share"] or 0 for m in ms)
        for m in sorted(ms, key=lambda m: -(m["share"] or 0)):
            share = m["share"] or 0.0
            renorm = 100 * share / total if total else 0.0
            abund = float(m["mpa_species"] or 0.0)
            rows.append([s, m["name"], m["species"] or "NA", m["mpa_name"] or "not found in MetaPhlAn",
                         float(m["completeness"] or 0), float(m["contamination"] or 0), share, renorm, abund])
            x.append(abund)
            y.append(renorm)
    write(os.path.join(a.outdir, "mag_vs_metaphlan.tsv"),
          ["sample", "mag", "tormes_species", "metaphlan_species", "completeness", "contamination",
           "pct_microbial_reads", "pct_of_kept_mag_reads", "metaphlan_pct"], rows)

    species = metaphlan["profiles"].get("species", [])
    cols = metaphlan["samples"]
    no_mag, dom = [], []
    for s in samples:
        i = cols.index(s)
        have = {(m["mpa_name"] or "").lower() for m in kept if m["sample"] == s}
        prof = sorted(((t["values"][i], t["taxon"]) for t in species if t["values"][i] > 0), reverse=True)
        for v, name in prof:
            if v >= a.min_abundance and name.lower() not in have:
                no_mag.append([s, name, float(v)])
        ms = sorted((m for m in kept if m["sample"] == s), key=lambda m: -(m["share"] or 0))
        top_mag = (ms[0]["mpa_name"] or ms[0]["species"] or "NA") if ms else "none"
        top_mpa = prof[0][1] if prof else "none"
        dom.append([s, top_mpa, float(prof[0][0]) if prof else 0.0, top_mag,
                    "yes" if top_mag.lower() == top_mpa.lower() else "NO"])
    write(os.path.join(a.outdir, "species_without_mag.tsv"), ["sample", "metaphlan_species", "metaphlan_pct"], no_mag)
    write(os.path.join(a.outdir, "dominant_species.tsv"),
          ["sample", "top_metaphlan", "metaphlan_pct", "top_kept_mag", "same"], dom)

    fate = []
    for s in samples:
        rows_d = mr.read_tsv(os.path.join(run, "10_metabat_depth", f"{s}.unmapped.depth.txt"))[1:]
        b = [(float(r[1]) * float(r[2]), float(r[1]), float(r[2])) for r in rows_d]
        tot = sum(v[0] for v in b) or 1.0
        hi = [v for v in b if v[2] > a.high_copy_depth]
        f = ig["fate"].get(s, {})
        fate.append([s, float(f.get("kept") or 0), float(f.get("other_bins") or 0), float(f.get("unbinned") or 0),
                     float(f.get("not_assembled") or 0), 100 * sum(v[0] for v in hi) / tot,
                     sum(v[1] for v in hi) / 1e6, len(hi)])
    write(os.path.join(a.outdir, "read_fate.tsv"),
          ["sample", "kept_mags", "other_bins", "unbinned", "not_assembled",
           "pct_mapped_on_high_copy_contigs", "high_copy_Mb", "high_copy_contigs"], fate)

    hum = []
    for s in samples:
        f = glob.glob(os.path.join(run, "04b_humann_output", f"{s}.*genefamilies.tsv"))
        if not f:
            continue
        un = tot = 0.0
        with open(f[0]) as h:
            for l in h:
                if l.startswith("#") or "|" in l.split("\t")[0]:
                    continue
                p = l.rstrip("\n").split("\t")
                tot += float(p[1])
                if p[0] == "UNMAPPED":
                    un = float(p[1])
        hum.append([s, un])
    write(os.path.join(a.outdir, "humann_unmapped.tsv"), ["sample", "humann_unmapped_reads"], hum)

    miss = [r for r in rows if r[8] == 0]
    lines = [f"run: {os.path.basename(run)}; samples: {len(samples)}",
             f"kept MAGs: {len(kept)}; with a species MetaPhlAn finds in the same sample: {len(kept) - len(miss)}",
             f"MAGs with no match: {', '.join(r[1] + ' (' + r[2] + ')' for r in miss) or 'none'}",
             f"MetaPhlAn species >= {a.min_abundance:g}% with no kept MAG: {len(no_mag)}",
             f"samples where the top species differs between branches: {sum(d[4] == 'NO' for d in dom)}/{len(dom)}"]
    if len(x) > 2:
        lines.append(f"MetaPhlAn vs share among kept MAGs, all kept MAGs: Pearson {pearson(x, y):.3f}, "
                     f"Spearman {pearson(ranks(x), ranks(y)):.3f}")
    with open(os.path.join(a.outdir, "summary.txt"), "w") as o:
        o.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

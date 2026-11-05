#!/usr/bin/env python3
"""Table of the screen: sample, group, read pairs, high-copy proxy, species, top species, and the sample
pair with the highest Bray-Curtis. Usage: screen_table.py <merged.txt> <kmer.tsv> <samplesheet.tsv>"""
import sys

rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1]) if not l.startswith("#")]
samples = [c.split(".")[0] for c in rows[0][1:]]
prof = {s: {} for s in samples}
for r in rows[1:]:
    if "s__" in r[0] and "t__" not in r[0]:
        for s, v in zip(samples, r[1:]):
            if float(v) > 0:
                prof[s][r[0].split("s__")[1]] = float(v)
kmer = {l.split()[0]: l.split() for l in open(sys.argv[2])}
sheet = {r[0]: r for r in (l.rstrip("\n").split("\t") for l in open(sys.argv[3]))}


def bray_curtis(a, b):
    return sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b)) / max(1e-9, sum(a.values()) + sum(b.values()))


print("sample\tgroup\tread_pairs\thigh_copy_pct\tspecies\tmost_different_sample\tbray_curtis\ttop_species")
for s in sorted(samples, key=lambda x: -float(kmer[x][2])):
    far = max((x for x in samples if x != s), key=lambda x: bray_curtis(prof[s], prof[x]))
    top = sorted(prof[s].items(), key=lambda x: -x[1])[:4]
    print(s, sheet[s][4], sheet[s][7], kmer[s][2], len(prof[s]), far, f"{bray_curtis(prof[s], prof[far]):.2f}",
          ", ".join(f"{n} {v:.0f}" for n, v in top), sep="\t")

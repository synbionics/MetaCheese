#!/usr/bin/env python3
"""Depth control: are a branch's blind spots real, or caused by low coverage?

Compares branch_tables.py results of the same samples at two depths. A species
that has a kept MAG at full depth but none at 5M pairs was a coverage effect.

Usage: compare_depth.py <tables_small_run> <tables_full_run>
"""
import os
import sys


def table(d, name):
    with open(os.path.join(d, name)) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]


small, full = sys.argv[1], sys.argv[2]
a, b = table(small, "species_without_mag.tsv"), table(full, "species_without_mag.tsv")
samples = {r["sample"] for r in table(small, "read_fate.tsv")} & {r["sample"] for r in table(full, "read_fate.tsv")}
gone = {(r["sample"], r["metaphlan_species"]) for r in b}
print("sample\tspecies\tmetaphlan_pct_small\tresult_at_full_depth")
for r in a:
    if r["sample"] in samples:
        key = (r["sample"], r["metaphlan_species"])
        print(r["sample"], r["metaphlan_species"], r["metaphlan_pct"],
              "still no kept MAG (not only coverage)" if key in gone else "kept MAG at full depth (coverage effect)", sep="\t")
print("\nsample\ttop_species_small_same\ttop_species_full_same", file=sys.stderr)
ds, df = {r["sample"]: r["same"] for r in table(small, "dominant_species.tsv")}, {r["sample"]: r["same"] for r in table(full, "dominant_species.tsv")}
for s in sorted(samples):
    print(s, ds.get(s), df.get(s), sep="\t", file=sys.stderr)

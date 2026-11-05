#!/usr/bin/env python3
"""Share of a sample's 31-mers that are high-copy, from the first reads of one FASTQ.

Prints: sample, then % of 31-mer occurrences with count >= 20, 40, 80 and 160, then
distinct 31-mers. One in 8 k-mers is counted, so run with PYTHONHASHSEED=0. On 500,000
reads, the >=40 column tracks the share of mapped bases on contigs above 2000x at 5M pairs.
"""
import collections
import subprocess
import sys

name, fq, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
p = subprocess.Popen(f"zcat '{fq}' | head -n {n * 4}", shell=True, stdout=subprocess.PIPE,
                     stderr=subprocess.DEVNULL, text=True)
c, k = collections.Counter(), 31
for i, l in enumerate(p.stdout):
    if i % 4 == 1:
        s = l.strip()
        for j in range(len(s) - k):
            m = s[j:j + k]
            if hash(m) & 7 == 0:
                c[m] += 1
tot = sum(c.values())
print(name, *[f"{sum(v for v in c.values() if v >= t) / tot * 100:.1f}" for t in (20, 40, 80, 160)], len(c), flush=True)

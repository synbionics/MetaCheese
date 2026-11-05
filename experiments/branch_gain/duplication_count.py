#!/usr/bin/env python3
"""Helper of duplication_source.sh: reads SAM lines on stdin, prints duplicate rates by contig depth.

Two measures, per class of contig depth:
- library: pairs with the same contig, position, strand and insert size. Two
  independent fragments rarely share all of these, so this mostly measures
  library (PCR, optical) duplication.
- by read: mate 1 reads whose first 50 bases were already seen in the whole
  file, close to what FastQC reports (it truncates long reads to 50 bases).
  Independent fragments of a very deep contig often start at the same base, so
  this also rises with real high-copy DNA.
"""
import sys

depth = {}
with open(sys.argv[1]) as f:
    next(f)
    for l in f:
        r = l.split("\t")
        depth[r[0]] = float(r[2])
CLASSES = [(0, 50), (50, 200), (200, 500), (500, 2000), (2000, float("inf"))]
n = len(CLASSES)
seen_pairs, total, dup, dup_read = [set() for _ in range(n)], [0] * n, [0] * n, [0] * n
seen_reads = set()
for l in sys.stdin:
    p = l.split("\t", 10)
    d = depth.get(p[2])
    if d is None:
        continue
    i = next(k for k, (a, b) in enumerate(CLASSES) if a <= d < b)
    rev = int(p[1]) & 16
    key = (p[2], p[3], rev, p[8])
    total[i] += 1
    if key in seen_pairs[i]:
        dup[i] += 1
    else:
        seen_pairs[i].add(key)
    start = (rev, p[9][-50:] if rev else p[9][:50])  # the read's own first 50 bases
    if start in seen_reads:
        dup_read[i] += 1
    else:
        seen_reads.add(start)
pct = lambda a, b: 100 * a / b if b else 0.0
print(f"{sys.argv[2]}: contig depth | pairs | % of pairs | library dup % | by-read dup %")
for k, (a, b) in enumerate(CLASSES):
    if total[k]:
        print(f"  {a:>5g}-{b:<6g} {total[k]:>9d} {pct(total[k], sum(total)):6.1f} {pct(dup[k], total[k]):6.1f} {pct(dup_read[k], total[k]):6.1f}")
print(f"  all          {sum(total):>9d}        {pct(sum(dup), sum(total)):6.1f} {pct(sum(dup_read), sum(total)):6.1f}")

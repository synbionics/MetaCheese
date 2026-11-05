#!/usr/bin/env python3
"""Add UniProt protein and organism names to a DIAMOND table against UniRef90. Records deleted from UniProt stay unnamed."""
import re
import sys
import urllib.request

rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
acc = lambda s: re.search(r"UniRef90_([A-Za-z0-9_-]+)", s).group(1)
ids = sorted({acc(r[1]) for r in rows if re.fullmatch(r"[A-Z0-9]{6,10}", acc(r[1]))})
names = {}
for i in range(0, len(ids), 40):
    q = "+OR+".join("accession:" + x for x in ids[i:i + 40])
    url = f"https://rest.uniprot.org/uniprotkb/search?fields=accession,protein_name,organism_name&format=tsv&size=50&query={q}"
    try:
        for l in urllib.request.urlopen(url, timeout=60).read().decode().splitlines()[1:]:
            a, p, o = l.split("\t")
            names[a] = (p, o)
    except OSError as e:
        print(f"UniProt request failed: {e}", file=sys.stderr)
print("protein\tuniref90\tpct_identity\tbitscore\tname\torganism")
for r in sorted(rows, key=lambda r: -float(r[7])):
    p, o = names.get(acc(r[1]), ("", ""))
    print(r[0], r[1].split("|")[0], r[2], r[7], p, o, sep="\t")

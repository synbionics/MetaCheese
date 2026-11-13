#!/usr/bin/env python3
"""Compare the summaries of two runs of the regression test.

Usage:
    python3 tests/regression/compare.py <baseline/summary.json> <new/summary.json>
        [--limits limits.json] [--out comparison.md] [--suggest-limits]

Each result is one of:
    FAIL   something broke: a failed run, a missing sample, step or file, or
           table columns that changed
    CHECK  the workflow ran, but a value moved beyond its limit, a species
           appeared or disappeared, or a version changed; a person decides
           whether the change makes sense
    PASS   within the limit

Exit code: 0 all pass, 1 at least one FAIL, 2 at least one CHECK.
With --suggest-limits, print the limits this pair of runs would need (use it on
two runs of the same version) and exit 0.
"""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Difference used for each per-sample value: "relative" is a fraction of the
# larger value, "absolute" is in the value's own unit.
PER_SAMPLE = {
    "trim_pairs": "relative", "retained_reads": "relative",
    "microbial_pairs": "relative", "host_pairs": "relative",
    "shannon": "absolute", "richness": "relative",
    "pathway_count": "relative", "gene_families": "relative",
    "humann_unmapped_pct": "absolute", "humann_unintegrated_pct": "absolute",
    "contigs": "relative", "assembly_bp": "relative", "n50": "relative",
    "mapping_rate": "absolute",
    "mags_total": "absolute", "mags_kept": "absolute",
    "mag_mean_completeness": "absolute", "mag_mean_contamination": "absolute",
    "kept_read_share": "absolute", "kept_mags_metaphlan_species": "absolute",
}
# Steps renamed since a stored baseline was made (old name -> new name): a step
# counts as completed under either name.
STEP_RENAMES = {"16.report.sh": "17.report.sh"}
# A species at this abundance (%) or above that is absent from the other run.
MAJOR_SPECIES_PCT = 1.0
# Used when turning observed noise into a limit, so that a value that happened
# not to move at all still gets some room where small changes are expected.
LIMIT_MARGIN = 2
EXACT_METRICS = {"trim_pairs", "retained_reads", "microbial_pairs", "host_pairs"}


def bray_curtis(a, b):
    """Bray-Curtis dissimilarity of two {feature: value} profiles (each scaled to sum 1)."""
    total_a, total_b = sum(a.values()), sum(b.values())
    if not total_a and not total_b:
        return 0.0
    if not total_a or not total_b:
        return 1.0
    return 0.5 * sum(abs(a.get(k, 0) / total_a - b.get(k, 0) / total_b) for k in set(a) | set(b))


def difference(kind, old, new):
    if kind == "absolute":
        return abs(old - new)
    largest = max(abs(old), abs(new))
    return abs(old - new) / largest if largest else 0.0


def compare(base, new, limits):
    """Return [{status, metric, item, baseline, new, change, limit}]."""
    rows = []

    def add(status, metric, item, old="", now="", change="", limit=""):
        rows.append({"status": status, "metric": metric, "item": item, "baseline": old,
                     "new": now, "change": change, "limit": limit})

    def value(metric, item, old, now, kind="absolute", diff=None):
        change = diff if diff is not None else difference(kind, old, now)
        limit = limits.get(metric)
        status = "PASS" if limit is None or change <= limit + 1e-12 else "CHECK"
        add(status, metric, item, old, now, round(change, 6), limit)

    if new.get("pipeline_exit") not in (0, None):
        add("FAIL", "pipeline", "run", "", f"exit {new['pipeline_exit']}")

    for s in base["samples"]:
        if s not in new["samples"]:
            add("FAIL", "sample missing", s)
    for s in new["samples"]:
        if s not in base["samples"]:
            add("CHECK", "new sample", s)
    for step in base.get("steps", []):
        if step not in new.get("steps", []) and STEP_RENAMES.get(step) not in new.get("steps", []):
            add("FAIL", "step not completed", step)
    for path, present in base["files"].items():
        if present and not new["files"].get(path):
            add("FAIL", "file missing or empty", path)
    for path, cols in base.get("columns", {}).items():
        if path in new.get("columns", {}) and new["columns"][path] != cols:
            add("FAIL", "table columns changed", path, ", ".join(cols), ", ".join(new["columns"][path]))

    for name in sorted(set(base.get("versions", {})) | set(new.get("versions", {}))):
        old, now = base.get("versions", {}).get(name), new.get("versions", {}).get(name)
        if old != now:
            add("CHECK", "version changed", name, old, now)

    shared = [s for s in base["samples"] if s in new["samples"]]
    for s in shared:
        a, b = base["per_sample"][s], new["per_sample"][s]
        for metric, kind in PER_SAMPLE.items():
            old, now = a.get(metric), b.get(metric)
            if old is None and now is None:
                continue
            if old is None or now is None:
                add("FAIL" if now is None else "CHECK", f"{metric} missing", s, old, now)
                continue
            value(metric, s, old, now, kind)
        value("species_bray_curtis", s, "", "", diff=bray_curtis(a["species"], b["species"]))
        value("pathway_bray_curtis", s, "", "", diff=bray_curtis(a["pathways"], b["pathways"]))
        for taxon, pct in sorted(a["species"].items()):
            if pct >= MAJOR_SPECIES_PCT and taxon not in b["species"]:
                add("CHECK", "species lost", f"{s}: {taxon}", pct, 0)
        for taxon, pct in sorted(b["species"].items()):
            if pct >= MAJOR_SPECIES_PCT and taxon not in a["species"]:
                add("CHECK", "species gained", f"{s}: {taxon}", 0, pct)
        old_names, new_names = list(a.get("mag_species", [])), list(b.get("mag_species", []))
        for name in list(old_names):
            if name in new_names:
                old_names.remove(name)
                new_names.remove(name)
        for name in old_names:
            add("CHECK", "MAG species lost", f"{s}: {name}")
        for name in new_names:
            add("CHECK", "MAG species gained", f"{s}: {name}")

    for pair, old in base.get("between_samples", {}).items():
        now = new.get("between_samples", {}).get(pair)
        if now is not None:
            value("between_sample_bray_curtis", pair, old, now)
    return rows


def suggest_limits(rows):
    """The limit each metric would need for these rows, with a margin."""
    largest = {}
    for r in rows:
        if isinstance(r["change"], (int, float)):
            largest[r["metric"]] = max(largest.get(r["metric"], 0.0), r["change"])
    out = {}
    for metric, seen in sorted(largest.items()):
        if metric in EXACT_METRICS and seen == 0:
            out[metric] = 0
        else:
            # two significant figures, rounded up
            limit = seen * LIMIT_MARGIN
            if limit:
                scale = 10 ** (1 - math.floor(math.log10(limit)))
                limit = math.ceil(limit * scale) / scale
            out[metric] = limit
    return out


def to_markdown(base, new, rows):
    counts = {k: sum(r["status"] == k for r in rows) for k in ("FAIL", "CHECK", "PASS")}
    lines = [f"# Regression comparison: {base['run']} vs {new['run']}", "",
             f"- Baseline: `{base['run']}`, commit {base.get('git_commit')}, {base.get('generated')}",
             f"- New: `{new['run']}`, commit {new.get('git_commit')}, {new.get('generated')}",
             f"- Result: **{counts['FAIL']} failed, {counts['CHECK']} to check, {counts['PASS']} passed**", ""]
    for status, title in (("FAIL", "Failed"), ("CHECK", "To check by hand")):
        chosen = [r for r in rows if r["status"] == status]
        if not chosen:
            continue
        lines += [f"## {title}", "", "| What | Item | Baseline | New | Change | Limit |", "|---|---|---|---|---|---|"]
        for r in chosen:
            lines.append(f"| {r['metric']} | {r['item']} | {r['baseline']} | {r['new']} | {r['change']} | {r['limit']} |")
        lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("baseline")
    ap.add_argument("new")
    ap.add_argument("--limits", default=os.path.join(HERE, "limits.json"))
    ap.add_argument("--out", help="also write the comparison to this Markdown file")
    ap.add_argument("--suggest-limits", action="store_true")
    args = ap.parse_args()

    with open(args.baseline) as fh:
        base = json.load(fh)
    with open(args.new) as fh:
        new = json.load(fh)
    with open(args.limits) as fh:
        limits = {k: v for k, v in json.load(fh).items() if not k.startswith("_")}

    rows = compare(base, new, limits)
    if args.suggest_limits:
        print(json.dumps(suggest_limits(rows), indent=1))
        return 0
    text = to_markdown(base, new, rows)
    print(text)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    if any(r["status"] == "FAIL" for r in rows):
        return 1
    return 2 if any(r["status"] == "CHECK" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())

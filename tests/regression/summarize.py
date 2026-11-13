#!/usr/bin/env python3
"""Summarise a run folder into a small test_summary/ folder for comparison.

Usage: python3 tests/regression/summarize.py output/<run> [--pipeline-exit N]

Writes output/<run>/test_summary/summary.json and copies a few small result
tables to test_summary/tables/. The report's own parsers are reused, so a
parser bug would show up in both. Python 3 standard library only.
"""
import argparse
import datetime as dt
import glob
import json
import os
import shutil
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "scripts", "report"))
import make_report as mr  # noqa: E402

# Must exist and not be empty, per sample ({s}) and per run.
PER_SAMPLE_FILES = [
    "00_fastqc/{s}_1_fastqc.zip", "00_fastqc/{s}_2_fastqc.zip",
    "01_AdapterRemoval/{s}_cleaned_L1.fq.gz", "01_AdapterRemoval/{s}_cleaned_L2.fq.gz",
    "01_AdapterRemoval/{s}.settings",
    "03_bowtie2_output/{s}.unmapped.fq.1.gz", "03_bowtie2_output/{s}.unmapped.fq.2.gz",
    "03_bowtie2_output/{s}.bowtie2.log",
    "04_metaphlan_output/{s}.unmapped.txt",
    "04b_humann_output/{s}.unmapped.fq_genefamilies.tsv",
    "04b_humann_output/{s}.unmapped.fq_pathabundance.tsv",
    "04b_humann_output/{s}.unmapped.fq_pathcoverage.tsv",
    "05_spades_output/{s}.unmapped/contigs.fasta",
    "05_spades_output/contigs/filtered/{s}.*",
    "07_Bowtie_Index/{s}.*.1.bt2*",
    "08_mapping_coverage/{s}.unmapped.sorted.bam",
    "10_metabat_depth/{s}.*depth.txt",
]
RUN_FILES = [
    "00_fastqc/multiqc_report.html",
    "04_metaphlan_output/merged_abundance_table.txt",
    "04_metaphlan_output/diversity/merged_abundance_table_bray-curtis.tsv",
    "04b_humann_output/merged_tables/humann_joined_pathabundance_cpm.tsv",
    "04b_humann_output/merged_tables/merged_table_ko.tsv",
    "04b_humann_output/merged_tables/merged_table_ec.tsv",
    "06b_quast/transposed_report.tsv",
    "13_checkm2/quality_report.tsv",
    "15_tormes_MAGs/taxonomic_identification/taxonomic-identification-kraken2.txt",
    "16_integration/read_fate.tsv",
    "16_integration/mags.tsv",
    "16_integration/mag_functions.tsv",
    "16_integration/mag_pathways.tsv",
    "16_integration/mag_ecs.tsv",
    "16_integration/ec_coverage.tsv",
    "report/report.html",
]
# The columns of these tables must not change.
COLUMN_TABLES = {
    "04_metaphlan_output/merged_abundance_table.txt": "#",
    "13_checkm2/quality_report.tsv": None,
    "15_tormes_MAGs/taxonomic_identification/taxonomic-identification-kraken2.txt": None,
    "04b_humann_output/merged_tables/merged_table_ko.tsv": None,
    "04b_humann_output/merged_tables/merged_table_ec.tsv": None,
    "16_integration/read_fate.tsv": None,
    "16_integration/mags.tsv": None,
    "16_integration/mag_functions.tsv": None,
    "16_integration/mag_pathways.tsv": None,
    "16_integration/mag_ecs.tsv": None,
    "16_integration/ec_coverage.tsv": None,
}
# Kept as they are, for reading and for summarising again later.
COPIED_TABLES = [
    "04_metaphlan_output/merged_abundance_table.txt",
    "04_metaphlan_output/diversity/merged_abundance_table_bray-curtis.tsv",
    "04_metaphlan_output/diversity/merged_abundance_table_shannon.tsv",
    "04_metaphlan_output/diversity/merged_abundance_table_richness.tsv",
    "04b_humann_output/merged_tables/unstratified/humann_joined_pathabundance_cpm.tsv",
    "06b_quast/transposed_report.tsv",
    "13_checkm2/quality_report.tsv",
    "15_tormes_MAGs/taxonomic_identification/taxonomic-identification-kraken2.txt",
    "16_integration/read_fate.tsv",
    "16_integration/mags.tsv",
    "16_integration/mag_functions.tsv",
    "16_integration/ec_coverage.tsv",
    "step_times.tsv",
    "versions.tsv",
]


def has_file(run, pattern):
    return any(os.path.getsize(p) > 0 for p in glob.glob(os.path.join(run, pattern)) if os.path.isfile(p))


def humann_sample(run, sample):
    """Pathway abundances and the unmapped/unintegrated shares of one sample."""
    out = {"pathways": {}, "unmapped_pct": None, "unintegrated_pct": None, "gene_families": None}
    path = os.path.join(run, "04b_humann_output", f"{sample}.unmapped.fq_pathabundance.tsv")
    if os.path.isfile(path):
        special = {}
        for row in mr.read_tsv(path, comment="#"):
            if "|" in row[0] or len(row) < 2:
                continue
            value = mr.to_float(row[1]) or 0.0
            if row[0] in ("UNMAPPED", "UNINTEGRATED"):
                special[row[0]] = value
            elif value > 0:
                out["pathways"][row[0]] = round(value, 4)
        total = sum(out["pathways"].values()) + sum(special.values())
        if total:
            out["unmapped_pct"] = round(100 * special.get("UNMAPPED", 0) / total, 4)
            out["unintegrated_pct"] = round(100 * special.get("UNINTEGRATED", 0) / total, 4)
    path = os.path.join(run, "04b_humann_output", f"{sample}.unmapped.fq_genefamilies.tsv")
    if os.path.isfile(path):
        out["gene_families"] = sum(1 for r in mr.read_tsv(path, comment="#") if "|" not in r[0] and r[0] != "UNMAPPED")
    return out


def diversity_column(run, metric):
    path = os.path.join(run, "04_metaphlan_output", "diversity", f"merged_abundance_table_{metric}.tsv")
    if not os.path.isfile(path):
        return {}
    return {mr.clean_sample(r[0]): mr.to_float(r[1]) for r in mr.read_tsv(path)[1:] if len(r) > 1}


def between_samples(run):
    path = os.path.join(run, "04_metaphlan_output", "diversity", "merged_abundance_table_bray-curtis.tsv")
    if not os.path.isfile(path):
        return {}
    rows = mr.read_tsv(path)
    names = [mr.clean_sample(c) for c in rows[0]]
    pairs = {}
    for r in rows[1:]:
        a = mr.clean_sample(r[0])
        for b, v in zip(names, r[1:]):
            if a < b:
                pairs[f"{a} / {b}"] = round(mr.to_float(v) or 0.0, 5)
    return pairs


def summarize(run, pipeline_exit=None):
    config_text = mr.read_text(os.path.join(run, "config.yml"), "") or ""
    sheet, _ = mr.collect_samplesheet(run, config_text, None)
    trimming = mr.collect_trimming(run)
    host = mr.collect_host(run)
    samples = [s["sample"] for s in sheet] or sorted(set(trimming) | set(host))

    metaphlan = mr.collect_metaphlan(run, samples)
    mapping = mr.collect_mapping(run, samples)
    mags = mr.collect_mags(run, samples)
    assembly = mr.collect_assembly(run)
    taxonomy = mr.collect_tormes_taxonomy(run)
    integration = mr.collect_joint(run) or {}  # R.1 tables
    fate = integration.get("fate") or {}
    joint_mags = {m["name"]: m for m in integration.get("mags") or []}
    functions, ec_cov = integration.get("functions") or [], integration.get("ec_coverage") or {}
    shannon, richness = diversity_column(run, "shannon"), diversity_column(run, "richness")

    species = {s: {} for s in samples}
    if metaphlan:
        for t in metaphlan["profiles"].get("species", []):
            for s, v in zip(metaphlan["samples"], t["values"]):
                if v and s in species:
                    species[s][t["taxon"]] = round(v, 5)

    per_sample = {}
    for s in samples:
        t, h, a = trimming.get(s, {}), host.get(s, {}), assembly.get(s, {})
        hum = humann_sample(run, s)
        kept = [m for m in mags if m["sample"] == s and m["passed"]]
        microbial, pairs_in = h.get("not_concordant"), h.get("pairs_in")
        mean = lambda key: round(statistics.mean(m[key] for m in kept), 3) if kept and all(m[key] is not None for m in kept) else None
        per_sample[s] = {
            "trim_pairs": t.get("pairs"),
            "retained_reads": t.get("retained_reads"),
            "microbial_pairs": microbial,
            "host_pairs": pairs_in - microbial if pairs_in is not None and microbial is not None else None,
            "shannon": shannon.get(s),
            "richness": richness.get(s),
            "species": species[s],
            "pathways": hum["pathways"],
            "pathway_count": len(hum["pathways"]),
            "gene_families": hum["gene_families"],
            "humann_unmapped_pct": hum["unmapped_pct"],
            "humann_unintegrated_pct": hum["unintegrated_pct"],
            "contigs": a.get("contigs"),
            "assembly_bp": a.get("total_bp"),
            "n50": a.get("n50"),
            "mapping_rate": mapping.get(s),
            "mags_total": sum(1 for m in mags if m["sample"] == s),
            "mags_kept": len(kept),
            "mag_mean_completeness": mean("completeness"),
            "mag_mean_contamination": mean("contamination"),
            "kept_read_share": round(fate[s]["kept"], 3) if s in fate and fate[s]["kept"] is not None else None,
            "kept_mags_metaphlan_species": sum(joint_mags.get(m["name"], {}).get("mpa_species") is not None for m in kept)
                                           if integration else None,
            "kept_mags_humann_species": sum(1 for r in functions if r["sample"] == s and r["humann_species"])
                                        if functions else None,
            "kept_mags_ec_share": round(ec_cov[s]["share"], 3) if ec_cov.get(s, {}).get("share") is not None else None,
            "mag_species": sorted((taxonomy.get(m["name"], {}).get("species") or "unclassified") for m in kept),
        }

    columns = {}
    for rel, comment in COLUMN_TABLES.items():
        path = os.path.join(run, rel)
        if os.path.isfile(path):
            rows = mr.read_tsv(path, comment=comment)
            if rows:
                columns[rel] = rows[0]

    versions = {}
    path = os.path.join(run, "versions.tsv")
    if os.path.isfile(path):
        versions = {r[0]: r[1] for r in mr.read_tsv(path) if len(r) > 1}

    times = {}
    path = os.path.join(run, "step_times.tsv")
    if os.path.isfile(path):
        times = {r[0]: mr.to_float(r[1]) for r in mr.read_tsv(path) if len(r) > 1}

    commit = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    return {
        "run": os.path.basename(os.path.normpath(run)),
        "generated": dt.date.today().isoformat(),
        "git_commit": commit,
        "pipeline_exit": pipeline_exit,
        "samples": samples,
        "files": {
            **{p.format(s=s): has_file(run, p.format(s=s)) for p in PER_SAMPLE_FILES for s in samples},
            **{p: has_file(run, p) for p in RUN_FILES},
        },
        "columns": columns,
        "steps": sorted(times),
        "step_seconds": times,
        "versions": versions,
        "per_sample": per_sample,
        "between_samples": between_samples(run),
        "mags": [{"name": m["name"], "sample": m["sample"], "kept": m["passed"],
                  "completeness": m["completeness"], "contamination": m["contamination"],
                  "size": m["size"], "species": taxonomy.get(m["name"], {}).get("species"),
                  "metaphlan_name": joint_mags.get(m["name"], {}).get("mpa_name"),
                  "metaphlan_species_pct": joint_mags.get(m["name"], {}).get("mpa_species")} for m in mags],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run folder, output/<run>")
    ap.add_argument("--pipeline-exit", type=int, default=None, help="exit code of the pipeline run")
    args = ap.parse_args()
    if not os.path.isdir(args.run):
        sys.exit(f"ERROR: {args.run} is not a folder")
    out = os.path.join(args.run, "test_summary")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(os.path.join(out, "tables"))
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summarize(args.run, args.pipeline_exit), fh, indent=1, sort_keys=True)
        fh.write("\n")
    for rel in COPIED_TABLES:
        path = os.path.join(args.run, rel)
        if os.path.isfile(path):
            shutil.copy(path, os.path.join(out, "tables", os.path.basename(rel)))
    print(f"Summary written to {out}")


if __name__ == "__main__":
    main()

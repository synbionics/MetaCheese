# Full example: PRJNA1198075, 28 cheeses

The complete dataset analysed in the MetaCheese manuscript: shotgun metagenomes
of 28 Grana-type hard cheeses, 13 PDO Grana Padano (`D_GP_*`) and 15 non-PDO (`ND_GP_*`),
Illumina NovaSeq 6000 paired-end 150 bp, 30 to 85 million read pairs per sample.
BioProject [PRJNA1198075](https://www.ebi.ac.uk/ena/browser/view/PRJNA1198075).

Try the [small example](../small_PRJNA1198075/) first: it uses the same
scripts on 6 samples and finishes in hours instead of days.

## Requirements

- the `metacheese` Docker image (see the main [README](../../README.md#quick-start));
- a server with **≥ 32 cores and ≥ 128 GB RAM** ([`Config.yml`](Config.yml) uses
  32 threads and 120 GB for metaSPAdes; adjust `par*` keys to your machine);
- about **1 TB of free disk**: 125 GB of FASTQ, plus intermediate files
  (trimmed and host-filtered reads, assemblies, BAM files);
- the *Bos taurus* genome in `data/gene/Bos_taurus/` (see the main README).

## Run it

From the repository root:

```bash
# 1. Download all FASTQ files from ENA, with metadata and sample groups
#    (resumes interrupted downloads and verifies every file's MD5)
bash examples/PRJNA1198075/download_PRJNA1198075.sh PRJNA1198075 /data/PRJNA1198075_data 4

# 2. Link (or copy) them into input/PRJNA1198075/ with the pipeline's naming
bash examples/PRJNA1198075/prepare_input.sh /data/PRJNA1198075_data

# 3. Run (output in output/example_PRJNA1198075/); use nohup/tmux, it takes days
nohup bash examples/PRJNA1198075/run_example.sh > example_full.out 2>&1 &
```

The download script writes `samplesheet.tsv` (sample name, run accession,
BioSample, group, file names, read pairs); the report uses its `group` column
to colour PDO and non-PDO samples.

If a step fails, fix the cause and resume with
`bash examples/PRJNA1198075/run_example.sh --from <step>` (e.g. `--from 05`).

## What to expect

- **Time**: HUMAnN (DIAMOND search) and metaSPAdes dominate. Scaling the small
  example's times by the read depth (×6 to ×17 per sample) and the number of
  samples (×4.7) gives an order of magnitude of several days on 32 cores.
- **Memory**: metaSPAdes is the only step that may need more than 64 GB on
  deep samples; it stops with an error when it reaches `step05.par2`.
- **MAGs**: with full depth, the dominant lactic acid bacteria (*Lactobacillus
  delbrueckii*, *L. helveticus*, *Lacticaseibacillus rhamnosus*, *Streptococcus
  thermophilus*) are expected to give medium- and high-quality MAGs.

## Notes on the data

- ENA serves `<run>_1.fastq.gz` / `<run>_2.fastq.gz`; `prepare_input.sh`
  renames them to `<sample>_1.fq.gz` / `<sample>_2.fq.gz`.
- The manuscript mentions 30 cheeses; the BioProject contains 28 runs.

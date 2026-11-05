# Small example: 6 cheeses, 5M read pairs each

A quick end-to-end test of MetaCheese on a workstation, built from the public
dataset [PRJNA1198075](https://www.ebi.ac.uk/ena/browser/view/PRJNA1198075)
(Grana-type hard cheeses: PDO Grana Padano and non-PDO). It runs every step of the workflow
and ends with the interactive report, so it is also the easiest way to see what
each step produces.

| Sample | Run | Group |
|--------|-----|-------|
| D_GP_02 | SRR31729973 | PDO |
| D_GP_04 | SRR31729955 | PDO |
| D_GP_13 | SRR31729948 | PDO |
| ND_GP_16 | SRR31729960 | non-PDO |
| ND_GP_10 | SRR31729968 | non-PDO |
| ND_GP_15 | SRR31729961 | non-PDO |

**How the samples were chosen.** By hand, not at random: the set covers both
very different and very similar communities, so that the report's diversity
sections have something to show. In this run, D_GP_04 / ND_GP_16 and
D_GP_02 / ND_GP_10 are well separated (Bray-Curtis 0.89 and 0.78), while
D_GP_13 / ND_GP_15 are similar (0.06). ND_GP_16 replaced ND_GP_4 after a first
run of this example, because ND_GP_4 fails FastQC's duplication check (60%).

Because the samples were selected for these contrasts, **this set cannot be
used to compare PDO with non-PDO cheeses**; use the
[full example](../PRJNA1198075/README.md) for that.

Only the **first 5 million read pairs** of each sample are kept (the full
samples have 30-85 million), so the profiles are noisier and MAGs are fewer and
less complete; the point is to exercise the pipeline, not to reproduce the
paper's results.

## Requirements

- the `metacheese` Docker image (see the main [README](../../README.md#quick-start));
- 16 cores and 64 GB RAM (lower `par*` threads and `step05.par2` in
  [`Config.yml`](Config.yml) for smaller machines);
- about 25 GB of free disk: 3.5 GB of input reads and 21 GB of results;
- the *Bos taurus* genome in `data/gene/Bos_taurus/` (see the main README, Quick start step 2;
  the index is built automatically on the first run).

## Run it

From the repository root, on the host:

```bash
# 1. Download the 6 samples, first 5M read pairs each -> input/small_PRJNA1198075/
bash examples/small_PRJNA1198075/prepare_input.sh

# 2. Run the whole pipeline (output in output/example_small/)
bash examples/small_PRJNA1198075/run_example.sh
```

`prepare_input.sh` streams each FASTQ file from ENA and stops after the first
5 million read pairs, so only about 3.5 GB are downloaded instead of the full
125 GB. Its speed depends on ENA: from a few minutes to a few hours. If the
connection drops it retries, and when run again it skips the files that are
already complete. Other options:

```bash
READ_PAIRS=1000000 bash examples/small_PRJNA1198075/prepare_input.sh   # a smaller test
bash examples/small_PRJNA1198075/prepare_input.sh /path/PRJNA1198075_data  # from the full download
```

The second form reads the files downloaded for the
[full example](../PRJNA1198075/README.md) and gives exactly the same reads.

`run_example.sh` passes extra options on to the runner:

```bash
bash examples/small_PRJNA1198075/run_example.sh --to 04.metaphlan    # stop after C.1/C.2
bash examples/small_PRJNA1198075/run_example.sh --from 05            # resume from the assembly
RUN_NAME=test2 bash examples/small_PRJNA1198075/run_example.sh       # a second, separate run
```

## What you get

`output/example_small/report/report.html` is the final report: open it in a
browser (it needs no internet connection). To send it to a colleague, use
`output/example_small/report/example_small_report.zip`, which also holds the
pages the report links to: unzip it and open `report.html`. The report explains
each step and shows its results; the other
folders are described in the main README ([Results](../../README.md#results)).

## Run time and resources

Measured on a 16-core / 62 GB workstation (Docker, local SSD/HDD):

| Step | Tool | Time |
|------|------|------|
| 00 | FastQC + MultiQC | 2 min |
| 00-01 | AdapterRemoval | 5 min |
| 03 | Bowtie2 host removal | 21 min |
| 04 | MetaPhlAn + diversity | 48 min |
| 04b-1 | HUMAnN (about 41 min per sample) | 4 h 8 min |
| 04b-last | HUMAnN join / renorm / regroup | < 1 min |
| 05 | metaSPAdes | 2 h 3 min |
| 06-07 | contig filter, QUAST, Bowtie2 index | 2 min |
| 08-09 | Bowtie2 mapping to contigs | 51 min |
| 10-11 | MetaBAT2 | < 1 min |
| 12-13 | CheckM2 | 5 min |
| 14 | MAG filter | < 1 min |
| 15 | TORMES | 34 min |
| 16 | report | < 1 min |
| **Total** | | **about 9 h** |

Peak memory was set by metaSPAdes (`step05.par2`) and HUMAnN (DIAMOND block
size `step04b_1.par2`). The run folder takes 21 GB.

## Expected results

All numbers come from `output/example_small/report/report.html`.

| Sample | Group | Microbial pairs (P.2) | Species (C.1) | Shannon | Top species (C.1) | Pathways (C.3) | Contig N50 | Kept MAGs | Reads in kept MAGs |
|--------|-------|------:|---:|---:|---|---:|---:|---:|---:|
| D_GP_02 | PDO | 4,869,435 | 14 | 1.30 | *L. delbrueckii* 59% | 142 | 9,095 | 4 | 21% |
| D_GP_04 | PDO | 4,741,276 | 28 | 1.21 | *L. delbrueckii* 60% | 320 | 2,123 | 4 | 58% |
| D_GP_13 | PDO | 4,949,888 | 15 | 0.67 | *Lacticaseibacillus rhamnosus* 86% | 146 | 2,921 | 4 | 68% |
| ND_GP_16 | non-PDO | 4,977,741 | 11 | 1.11 | *Lacticaseibacillus zeae* 65% | 210 | 10,005 | 6 | 76% |
| ND_GP_10 | non-PDO | 4,674,647 | 29 | 1.70 | *L. helveticus* 45% | 183 | 13,054 | 6 | 57% |
| ND_GP_15 | non-PDO | 4,498,918 | 19 | 0.60 | *Lacticaseibacillus rhamnosus* 87% | 144 | 5,545 | 5 | 65% |

- Host (*Bos taurus*) reads are 0.5-11% of each sample.
- HUMAnN maps 82-86% of the reads to UniRef90 genes in every sample. 324
  MetaCyc pathways are found across the six samples.
- MetaBAT2 makes 45 bins; 29 pass CheckM2 (completeness > 50%,
  contamination < 10%) and are annotated by TORMES. 26 of the 29 belong to a
  species that MetaPhlAn also finds in the same sample: *Lactobacillus*
  *delbrueckii* and *helveticus*, *Lacticaseibacillus rhamnosus* and
  *paracasei*, *Lentilactobacillus buchneri* and *parabuchneri*,
  *Schleiferilactobacillus harbinensis*, *Streptococcus thermophilus*.
- The other 3 match MetaPhlAn only at genus level: TORMES calls them
  *L. casei* (ND_GP_16, 63% of the sample's reads), *L. kefiri* (ND_GP_10, 16%)
  and *L. buchneri* (ND_GP_10, 2%). MetaPhlAn reports *Lacticaseibacillus zeae*
  65% and *Lentilactobacillus rapi* 21% in these samples. TORMES uses an older
  Kraken2 database, so these may be naming differences rather than real
  disagreements.
- In D_GP_02, kept MAGs hold only 21% of the reads, and the top species
  *L. delbrueckii* has no kept MAG: two bins holding 35% of the reads fail on
  completeness (48% and 37%) at this depth.
- TORMES' own summary report (`15_tormes_MAGs/*.html`) is not written: its R
  script fails on the MLST table when MAGs of species with and without an
  MLST scheme are mixed. The per-MAG results (annotation, MLST, resistance
  and virulence genes, `genome_stats/*/report.html`) are complete.

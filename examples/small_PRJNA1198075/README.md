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
| ND_GP_4 | SRR31729957 | non-PDO |
| ND_GP_10 | SRR31729968 | non-PDO |
| ND_GP_15 | SRR31729961 | non-PDO |

The six samples are those used as an example in the MetaCheese manuscript:
D_GP_04 / ND_GP_4 and D_GP_02 / ND_GP_10 are well separated in beta diversity,
while D_GP_13 and ND_GP_15 are similar. Only the **first 5 million read pairs**
of each sample are kept (the full samples have 30-85 million), so the profiles
are noisier and few MAGs reach good quality; the point is to exercise the
pipeline, not to reproduce the paper's results.

## Requirements

- the `metacheese` Docker image (see the main [README](../../README.md#quick-start));
- 16 cores and 64 GB RAM (lower `par*` threads and `step05.par2` in
  [`Config.yml`](Config.yml) for smaller machines);
- about 25 GB of free disk: 3.5 GB of input reads and 19 GB of results;
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
browser, or send the file to a colleague (it is a single file and needs no
internet connection). It explains each step and shows its results; the other
folders are described in the main README ([Results](../../README.md#results)).

## Run time and resources

Measured on a 16-core / 62 GB workstation (Docker, local SSD/HDD):

| Step | Tool | Time |
|------|------|------|
| 00 | FastQC + MultiQC | 1 min |
| 00-01 | AdapterRemoval | 20 min |
| 03 | Bowtie2 host removal | 25 min |
| 04 | MetaPhlAn + diversity | 59 min |
| 04b-1 | HUMAnN (about 37 min per sample) | 3 h 40 min |
| 04b-last | HUMAnN join / renorm / regroup | < 1 min |
| 05 | metaSPAdes | 1 h 57 min |
| 06-07 | contig filter, Bowtie2 index | 1 min |
| 08-09 | Bowtie2 mapping to contigs | 50 min |
| 10-11 | MetaBAT2 | < 1 min |
| 12-13 | CheckM2 | 5 min |
| 14 | MAG filter | < 1 min |
| 15 | TORMES | 32 min |
| 16 | report | < 1 min |
| **Total** | | **about 9 h** |

The run was resumed twice inside HUMAnN, so the HUMAnN time is extrapolated
from the 5 samples timed in one go; MetaPhlAn ran with 4 threads (now 16 in
the config). Peak memory was set by metaSPAdes (`step05.par2`) and HUMAnN
(DIAMOND block size `step04b_1.par2`). The run folder takes 19 GB.

## Expected results

All numbers come from `output/example_small/report/report.html`.

| Sample | Group | Microbial pairs (P.2) | Species (C.1) | Shannon | Top species (C.1) | Pathways (C.3) | Contig N50 | Kept MAGs | Reads in kept MAGs |
|--------|-------|------:|---:|---:|---|---:|---:|---:|---:|
| D_GP_02 | PDO | 4,869,435 | 14 | 1.30 | *L. delbrueckii* 59% | 142 | 8,818 | 4 | 23% |
| D_GP_04 | PDO | 4,741,276 | 28 | 1.21 | *L. delbrueckii* 60% | 320 | 1,939 | 4 | 57% |
| D_GP_13 | PDO | 4,949,888 | 15 | 0.67 | *Lacticaseibacillus rhamnosus* 86% | 146 | 2,791 | 4 | 68% |
| ND_GP_4 | non-PDO | 4,752,018 | 9 | 0.86 | *L. helveticus* 79% | 180 | 5,382 | 4 | 19% |
| ND_GP_10 | non-PDO | 4,674,647 | 29 | 1.70 | *L. helveticus* 45% | 183 | 12,403 | 6 | 56% |
| ND_GP_15 | non-PDO | 4,498,918 | 19 | 0.60 | *Lacticaseibacillus rhamnosus* 87% | 144 | 5,245 | 5 | 65% |

- Host (*Bos taurus*) reads are 1-11% of each sample.
- HUMAnN maps 82-86% of the reads to UniRef90 genes in every sample except
  ND_GP_4 (52%). 323 MetaCyc pathways are found across the six samples.
- MetaBAT2 makes 41 bins; 27 pass CheckM2 (completeness > 50%,
  contamination < 10%) and are annotated by TORMES. All 27 belong to a
  species that MetaPhlAn also finds in the same sample: *Lactobacillus*
  *delbrueckii*, *L. helveticus*, *Lacticaseibacillus rhamnosus* and
  *paracasei*, *Lentilactobacillus buchneri*, *parabuchneri* and *kefiri*,
  *Schleiferilactobacillus harbinensis*, *Streptococcus thermophilus*,
  *Anoxybacillus flavithermus*.
- TORMES' own summary report (`15_tormes_MAGs/*.html`) is not written: its R
  script fails on the MLST table when MAGs of species with and without an
  MLST scheme are mixed. The per-MAG results (annotation, MLST, resistance
  and virulence genes, `genome_stats/*/report.html`) are complete.

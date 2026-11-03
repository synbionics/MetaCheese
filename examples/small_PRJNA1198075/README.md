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
- about 50 GB of free disk for inputs (3.5 GB) and outputs;
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

<!-- RUNTIME_TABLE -->

## Expected results

<!-- RESULTS_SUMMARY -->

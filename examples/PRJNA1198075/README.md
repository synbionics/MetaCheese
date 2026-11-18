# Full example: PRJNA1198075, 28 cheeses

The complete dataset analysed in the MetaCheese manuscript: shotgun metagenomes
of 28 Grana-type hard cheeses, 13 PDO Grana Padano (`D_GP_*`) and 15 non-PDO (`ND_GP_*`),
Illumina NovaSeq 6000 paired-end 150 bp, 30 to 85 million read pairs per sample.
BioProject [PRJNA1198075](https://www.ebi.ac.uk/ena/browser/view/PRJNA1198075).

Try the [small example](../small_PRJNA1198075/) first: it uses the same
scripts on 6 samples and finishes in hours instead of days.

## Requirements

- the `metacheese` Docker image, built from the repository's `Dockerfile`;
- a server with **≥ 32 cores and ≥ 128 GB RAM** ([`Config.yml`](Config.yml) uses
  32 threads and 120 GB for metaSPAdes; adjust `par*` keys to your machine);
- about **1 TB of free disk**: 125 GB of FASTQ, plus intermediate files
  (trimmed and host-filtered reads, assemblies, BAM files);
- the *Bos taurus* genome in `data/gene/Bos_taurus/` (`download_host_genome.sh`).

## Preparing another machine

Everything is rebuilt on the server from the repository: the image from the
`Dockerfile`, the reads and the host genome from their public archives. Nothing
has to be copied from another machine.

| What | Size | How |
|------|------|-----|
| Repository | small | `git clone` of the commit to run |
| Docker image `metacheese` | about 130 GB, several hours to build (mostly database downloads) | `docker build -t metacheese .` |
| Reads, `PRJNA1198075_data/` | 116 GB, 56 FASTQ files | `download_PRJNA1198075.sh` (checks every file's MD5) |
| Host genome, `data/gene/Bos_taurus/` | 2.8 GB (plus 3.5 GB of index, built on the first run in about 40 min) | `download_host_genome.sh` (checks NCBI's MD5 and the FASTA's) |

The build gives the same tools and databases as the image that passed the
regression test (`tests/`): every conda environment is installed from the
exact package list of that image (`docker/envs/`), the databases are
downloaded at fixed versions, and TORMES's extra files at fixed commits
(`docker/tormes/setup.sh`). `check_server.sh` compares the built image's tool
and database versions with the tested ones.

On the server, from a directory with room for the data (`/data` below):

```bash
git clone <repository URL> MetaCheese && cd MetaCheese
docker build -t metacheese . 2>&1 | tee build.log
bash examples/PRJNA1198075/download_host_genome.sh
bash examples/PRJNA1198075/download_PRJNA1198075.sh PRJNA1198075 /data/PRJNA1198075_data 4
bash examples/PRJNA1198075/prepare_input.sh /data/PRJNA1198075_data
```

Then:

1. Set threads and memory in [`Config.yml`](Config.yml) for the server: every
   `par*` thread count at most the number of cores, and metaSPAdes memory
   (`step05.par2`, GB) below the RAM, with room for the system. The deepest
   sample (ND_GP_2, 85M read pairs) needs the most.
2. `bash examples/PRJNA1198075/check_server.sh`: checks the image (tool and
   database versions against the tested ones), cores and memory against the
   config, free disk, the 28 samples and the host genome. It changes nothing.
   Fix every FAIL before starting; a WARN is for you to judge.

## Run it

From the repository root:

```bash
# After "Preparing another machine" and a check_server.sh without FAIL.
# Output in output/example_PRJNA1198075/; use nohup/tmux, it takes days.
nohup bash examples/PRJNA1198075/run_example.sh > example_full.out 2>&1 &
```

The download script writes `samplesheet.tsv` (sample name, run accession,
BioSample, group, file names, read pairs); the report uses its `group` column
to colour PDO and non-PDO samples.

If a step fails, fix the cause and resume with
`bash examples/PRJNA1198075/run_example.sh --from <step>` (e.g. `--from 05`).

## What to expect

- **Time**: HUMAnN (DIAMOND search) and metaSPAdes dominate. The 28 samples
  hold 1,030 million read pairs (30 to 85 million each), 34 times the small
  example, which took about 9 h on 16 cores. Scaled linearly, that is about
  300 h on 16 cores, or roughly a week on 32; metaSPAdes may scale worse.
- **Memory**: metaSPAdes is the only step that may need more than 64 GB on
  deep samples; it stops with an error when it reaches `step05.par2`.
- **Resuming**: `--from <step>` runs that step again for every sample. Steps
  00-01 (trimming) and 04b-1 (HUMAnN) skip the samples they have already
  finished; the others, including metaSPAdes (05), start again from the first
  sample.
- **MAGs**: with full depth, the dominant lactic acid bacteria (*Lactobacillus
  delbrueckii*, *L. helveticus*, *Lacticaseibacillus rhamnosus*, *Streptococcus
  thermophilus*) are expected to give medium- and high-quality MAGs.

## Notes on the data

- ENA serves `<run>_1.fastq.gz` / `<run>_2.fastq.gz`; `prepare_input.sh`
  renames them to `<sample>_1.fq.gz` / `<sample>_2.fq.gz`.
- The manuscript mentions 30 cheeses; the BioProject contains 28 runs.

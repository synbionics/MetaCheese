# Configuration and advanced use

## The config file

`config/Config.yml` has one section per step (`step00`, `step00_01`,
`step03`, ...). Inside a section:

- `dir*` / `var*` keys are **paths**, relative to `scripts/pipeline/` inside
  the container, where the step scripts run;
- `par*` keys are **parameters**: threads, memory, thresholds.

`@input_sample@` and `@output_folder@` in the paths are replaced with the
input folder and the run name when a run starts. The resolved config is saved
as `output/<run>/config.yml` and reused when the run is resumed, so later
edits to `config/Config.yml` do not change a run that already exists.

The examples ship their own configs, sized for their machines
(`examples/*/Config.yml`).

`config/configurator.html` is a form that writes this file: open it in a browser
(no internet or installation needed). Threads follow the cores you enter, and the
metaSPAdes memory limit is your RAM minus 8 GB; each value can be edited.

## Parameters

| Section | Key | Meaning | Default |
|---------|-----|---------|---------|
| `step00` | `par1` | FastQC threads (one file per thread) | 8 |
| `step00_01` | `par1`-`par6` | AdapterRemoval threads, min length, min quality, max Ns, trim 5', trim 3' | 16, 50, 30, 10, 2, 2 |
| `step03` | `dir3`, `par1` | host Bowtie2 index prefix, threads | `Bos_taurus`, 16 |
| `step04` | `par1` | MetaPhlAn threads | 16 |
| `step04b_1` | `par1`, `par2` | HUMAnN threads, DIAMOND block size (~6 GB RAM per unit) | 16, 2 |
| `step05` | `par1`, `par2` | metaSPAdes threads, memory limit (GB) | 16, 56 |
| `step06` | `filter1` | parallel jobs for contig filtering | 32 |
| `step06b` | `par1`, `par2` | QUAST threads, minimum contig length counted (bp) | 16, 500 |
| `step07` | `par1` | Bowtie2-build threads (contig indexes) | 16 |
| `step08_09` | `par1` | Bowtie2/samtools threads | 16 |
| `step10_11` | `par1`, `par2` | MetaBAT2 min contig length, threads | 1500, 16 |
| `step12_13` | `par2`, `par3` | pplacer threads, CheckM2 threads | 4, 16 |
| `step14` | `var5`, `var6` | MAG filter: completeness above, contamination below (%) | 50, 10 |
| `step15` | `par2` | TORMES threads | 16 |

The meaning of every key of every step is in [`steps.md`](steps.md).

## Speed

Every multi-threaded tool takes its thread count from the config, so the main
lever is matching the `par*` values to your cores. The pipeline already uses
options that are faster without changing the results (HUMAnN
`--memory-use maximum`, AdapterRemoval `--gzip-level 1`). With spare RAM,
raise `step04b_1.par2` (DIAMOND block size) to speed up HUMAnN, the slowest
community-level step.

Two further options are faster but change the results slightly, so they are
**not** enabled:

- Bowtie2 `--very-fast` for host removal (template `03`): 2-3× faster, a few
  more host reads pass through.
- The EC-filtered UniRef90 database for HUMAnN
  (`humann_databases --download uniref uniref90_ec_filtered_diamond`): 5-10×
  faster translated search, but only enzymes are found, so part of the gene
  families and KO groups are lost.

## How templates and config are connected

`scripts/run_pipeline.sh` (function `generate_scripts`) copies each
`scripts/templates/template_<step>.sh` to `scripts/pipeline/` and replaces the
placeholders `@<prefix>_<name>@` with the value of `<step section>.<key>` from
the config. The mapping is the `STEPS` table at the top of `run_pipeline.sh`:

```
[step04b_1]="template_04b-1.humann_job_array.sh  scripts/pipeline/04b-1.humann_job_array.sh  04b_1  var1=dir1 ... par2=par2"
              template                           generated script                            prefix placeholder=config key
```

so `@04b_1_par2@` in the template becomes the value of `step04b_1.par2`.

To add a parameter: use a new placeholder in the template, add `name=key` to
the step's `STEPS` entry, and add the key to every config
(`config/Config.yml`, `examples/*/Config.yml`). Also update the copy of
`config/Config.yml` inside `config/configurator.html` (block `id="template"`), and
add a field to its `GROUPS` table if users should set the key. A key missing from the config
stops the run with an error.

The scripts in `scripts/pipeline/` are regenerated at the start of each run:
do not regenerate them while a run is using them.

## Other ways to run

### Inside a container you started yourself

```bash
docker compose up -d                       # uses docker-compose.yml
docker exec -it metacheese-container bash
```

or

```bash
docker run -it --rm -v "$PWD/scripts:/main/scripts" -v "$PWD/config:/main/config" \
    -v "$PWD/data:/main/data" -v "$PWD/input:/main/input" -v "$PWD/output:/main/output" \
    -w /main metacheese bash
```

then, inside the container:

```bash
bash scripts/run_pipeline.sh --input my_dataset --name my_run --config config/Config.yml
```

with the same `--from` / `--to` options as `examples/run_in_docker.sh`.
`scripts/`, `config/`, `data/`, `input/` and `output/` are mounted from the
repository, so edits to templates and configs take effect without rebuilding
the image. Files written by the container belong to root;
`examples/run_in_docker.sh` gives them back to your user at the end, here you
have to do it yourself (`sudo chown -R $USER output/<run>`).

### Interactive menu

`bash scripts/run_pipeline.sh` without arguments (inside the container) opens
a menu: start a new run (it asks for the input folder and a code, and creates
`output/<date>_<code>/`) or continue an existing one from a given step.

### Regenerating only the report

```bash
python3 scripts/report/make_report.py output/<run>     # on the host, Python 3 only
```

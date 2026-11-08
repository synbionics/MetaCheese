# Tests

Two kinds of test, for two different questions.

## Unit tests: is the code correct?

They need no Docker and no data, and take less than a second:

```bash
python3 -m unittest discover -s tests/unit
```

They cover the host checks of `examples/run_in_docker.sh` (using a stand-in for
Docker), the configs, the configurator's copy of `config/Config.yml`, the
report's parsers and the regression test's own code.

## Regression test: does the workflow still give sensible results?

It runs every step on samples 1, 2 and 3 of each group of
[PRJNA1198075](https://www.ebi.ac.uk/ena/browser/view/PRJNA1198075) (`D_GP_01`,
`D_GP_02`, `D_GP_03`, `ND_GP_1`, `ND_GP_2`, `ND_GP_3`), first 5M read pairs of
each, and compares the results with those stored in `regression/baseline/`. It
takes about 9 hours on 16 cores / 62 GB and uses about 21 GB in `output/`.

Run it after a change to a tool, a database or a step, one change at a time:

```bash
bash tests/regression/prepare_input.sh [FULL_DATA_DIR]   # once; 3.4 GB
bash tests/regression/run.sh my_update
```

Without `FULL_DATA_DIR` the reads are streamed from ENA. To test another image,
set `IMAGE=metacheese:some_tag`. To resume, add `--from <step>` as for any run.

The comparison is printed and saved in `output/<run>/test_comparison.md`:

| Result | Meaning | Exit code |
|--------|---------|-----------|
| **FAIL** | Something broke: a failed run, a missing sample, step or file, or table columns that changed | 1 |
| **CHECK** | It ran, but a value moved beyond its limit, a species appeared or disappeared, or a version changed. A person decides whether the change makes sense | 2 |
| **PASS** | Everything is within its limit | 0 |

An update to a database or tool is expected to give CHECK results, so CHECK does
not mean wrong. It lists what to look at.

### What is compared

- Per sample: read counts, Shannon and richness, species profile, HUMAnN pathways
  and unmapped shares, assembly size and N50, mapping rate, MAG counts and
  quality, share of reads in kept MAGs, species names of the kept MAGs.
- Between samples: Bray-Curtis.
- Per run: expected files, table columns, completed steps, tool and database
  versions.

MAGs are compared by count, quality and species name per sample, because bin
numbers are not stable between runs.

### Updating the stored results

`regression/baseline/` holds `summary.json` and a few small tables, not the run
itself. After reviewing and accepting a change, make that run the new baseline:

```bash
rm -rf tests/regression/baseline
cp -r output/<run>/test_summary tests/regression/baseline
```

`regression/limits.json` says how far each value may move. Set the limits from
two runs of the same version, which show how much results vary with no change at
all:

```bash
python3 tests/regression/compare.py tests/regression/baseline/summary.json \
    output/<repeat_run>/test_summary/summary.json --suggest-limits
```

The current image was tagged `metacheese:baseline-2026-11` so that the baseline
can be reproduced after a rebuild. Keep that tag until the new baseline is
accepted.

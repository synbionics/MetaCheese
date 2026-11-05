# Experiment: what does running both branches gain?

**Status: proposal, for review. Not committed, and not part of the pipeline.**
Numbers in this file are preliminary: 5M read pairs, from the archived run
`output/example_small_with_ND_GP_4`, with tables in `preliminary_results/`.

## Question

Community-level profiling (MetaPhlAn, HUMAnN) and assembly-specific profiling
(metaSPAdes, MetaBAT2, CheckM2, TORMES) are usually run separately. On the same
reads, what would a researcher conclude wrongly, or fail to see, using only one?

The experiment compares three readings of the same run. No new pipeline run
is needed for the comparison itself:

| Arm | Uses | Stands for |
|---|---|---|
| Community only | C.1-C.4 | A profiling study |
| Assembly only | A.1-A.7 | A MAG study |
| Joint | R.1 and these scripts | MetaCheese |

## Data

The 28-sample full-depth run (`examples/PRJNA1198075`), to be run on the larger
machine, plus the 6-sample small run as the low-depth comparison.

## Measures

| # | Measure | Script | Output | Preliminary (5M pairs, 6 samples) |
|---|---|---|---|---|
| M1 | Assembly-only blind spots: species >= 1% with no kept MAG; samples where the top species changes | `branch_tables.py` | `species_without_mag.tsv`, `dominant_species.tsv` | 13 species; top species differs in 1 of 6 samples (D_GP_02) |
| M2 | Community-only blind spots: reads on high-copy contigs, HUMAnN unmapped reads | `branch_tables.py`, `phage_check.sh` | `read_fate.tsv`, `humann_unmapped.tsv`, `named_hits.tsv` | ND_GP_4: 56% of mapped bases on 0.17 Mb; others 1.5-13% |
| M3 | Species names: TORMES vs MetaPhlAn, settled by ANI | `branch_tables.py`, `species_check.sh` | `mag_vs_metaphlan.tsv`, `bins_vs_references.tsv` | 2 of 27 kept MAGs disagree; Mash ANI supports MetaPhlAn in both |
| M4 | Agreement where both see a species: MetaPhlAn vs share among kept MAGs | `branch_tables.py` | `summary.txt` | Pearson 0.94, Spearman 0.78 (27 MAGs, including the 2 misnamed) |
| C1 | Depth control: is a blind spot real, or low coverage? | `compare_depth.py` | Printed table | Needs the full run |
| C2 | Duplication: real high-copy DNA, or library artefact? | `duplication_source.sh` | `duplication_by_depth.txt` | See below |

## The three cases so far

Each is one sample at 5M pairs, so these are leads to test at full depth.

- **ND_GP_4, community only.** MetaPhlAn assigns 96% of the profile to four
  species. 56% of mapped bases are on 0.17 Mb of contigs above 2000x. Several
  contigs carry phage genes, with closest hits at 30-60% identity, and the top contig
  (21% of mapped bases) is still unidentified. HUMAnN leaves about 3.1M reads
  unmapped here and 1.3-1.6M elsewhere. The link between those unmapped reads and
  these contigs is an inference, not yet measured. The bacterial proportions are not
  distorted: *L. helveticus* is 78.4% of kept-MAG reads and 78.6% in MetaPhlAn.
- **D_GP_02, assembly only.** MetaPhlAn gives *L. delbrueckii* 59%, with no kept
  MAG. Two rejected bins (48% and 37% complete, about 35% of reads) are closest
  to *L. delbrueckii*. Kept MAGs make *L. rhamnosus* look dominant. A rejected,
  mixed casei-group bin (99% complete, 38% contamination) fits MetaPhlAn's
  *L. paracasei* and *L. zeae*. May be a coverage effect: C1 decides.
- **ND_GP_10, naming.** TORMES' "*L. kefiri*" is closest to *L. rapi* (Mash
  94.5%, MetaPhlAn 21.5%), and "*L. buchneri*" to *L. otakiensis* (98.1%, MetaPhlAn
  5.4%). The likely cause is TORMES' 2020 minikraken database, so this is a
  database limit, not biological complementarity.

## Duplication in ND_GP_4

FastQC reports 60% duplication (FAIL) in ND_GP_4 and 20-33% in the others.
`duplication_by_depth.txt` separates two parts:

- **Library duplication is about 30-37% at every contig depth.** It is flat, so
  it does not come from high-copy DNA. It is similar to D_GP_13 (34%), and higher
  than the other four samples (14-23%).
- **The excess is on the deep contigs.** 56% of pairs are on contigs above
  2000x, where 84% of reads are exact duplicates by their first 50 bases.

So the FastQC failure combines a fairly high library rate with real high-copy DNA.
The by-read measure runs 6-12 points above FastQC (66% vs 60% for ND_GP_4), and
ranks the samples almost the same way. It uses trimmed, mapped reads.

## Scripts

Run from the repository root.

| Script | Purpose | Tested |
|---|---|---|
| `branch_tables.py <run> <outdir>` | M1-M4 tables | Yes, on the archived run |
| `duplication_source.sh <run> <samples...>` | C2 | Yes, on the 6 samples |
| `species_check.sh <run> <outdir>` | Bins vs reference genomes | Yes, on the archived run |
| `compare_depth.py <small> <full>` | C1 | Only logic: no full run yet |
| `phage_check.sh <run> <sample> <outdir>` | Annotates high-copy contigs | Commands ran by hand on ND_GP_4; the script itself has not run |
| `screen_samples.sh <data> <outdir>` | Quick screen of all samples | Commands ran by hand on 28 samples; the script itself has not run |
| `name_hits.py`, `kmer_highcopy.py`, `screen_table.py`, `duplication_count.py` | Helpers | Yes |

Typical use after the full run:

```bash
python3 experiments/branch_gain/branch_tables.py output/example_PRJNA1198075 results/full
python3 experiments/branch_gain/branch_tables.py output/example_small results/small
python3 experiments/branch_gain/compare_depth.py results/small results/full
bash experiments/branch_gain/species_check.sh output/example_PRJNA1198075 results/species
bash experiments/branch_gain/duplication_source.sh output/example_PRJNA1198075 ND_GP_4 ND_GP_16
bash experiments/branch_gain/phage_check.sh output/example_PRJNA1198075 ND_GP_4 results/phage
```

## Decisions for you

1. **Adopt this design?** It adds no heavy compute beyond the planned full run.
2. **Add tools to confirm the claims?** fastANI (species) and geNomad or CheckV
   (phage) are not in the image. Without them, these remain estimates and hints.
   They can go in a separate environment rather than the main image.
3. **Add a simulated community?** Known genomes, a phage and close casei-group
   species would give true precision and recall, which real data cannot give. It
   needs a simulator such as InSilicoSeq. Not scripted yet.
4. **Test shared preprocessing?** The README's claim has little support so far:
   Bray-Curtis between a sample's 1M raw and 5M filtered profiles is 0.01-0.02
   (MetaPhlAn only, not assemblies). Either test it with a second preprocessing
   run, or drop it from the main argument. Not scripted yet.
5. **Test deduplication?** The pipeline does not remove duplicates. Library rates
   range from 14% to 35%, which could matter for MAG depth and abundance.
   Not scripted yet.

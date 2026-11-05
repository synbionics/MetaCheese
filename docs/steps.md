# MetaCheese step reference

One section per step, in execution order. IDs refer to the workflow figure
([`img/workflow.svg`](img/workflow.svg)). `<s>` is the sample name (the input
file name without `_1.fq.gz`), `<run>` the run folder `output/<run>/`.
Config keys are under the step's section in the config file; paths are relative
to `scripts/pipeline/` inside the container.

Execution order (`STEPS_ORDERED` in `scripts/run_pipeline.sh`):

```
00 → 00-01 → 03 → 04 → 04b-1 → 04b-last → 05 → 06 → 06b → 07 → 08-09 → 10-11 → 12-13 → 14 → 15 → 16
```

---

## I.2: Host genome index

| | |
|---|---|
| Script | `scripts/build_bowtie2_index.sh [HOST] [THREADS]` (run once, not part of `STEPS_ORDERED`) |
| Tool | `bowtie2-build` (Bowtie2 2.3.3.1, `bowtieenv`) |
| Input | `data/gene/<HOST>/*.fa|*.fna|*.fasta` (downloaded by the user, see README) |
| Output | `data/gene/<HOST>/<HOST>.{1-4,rev.1,rev.2}.bt2` |
| Notes | Without arguments it asks for the host name. An existing index is kept. *Bos taurus* ARS-UCD2.0 takes ~40 min with 16 threads and ~4 GB of index files. |

## Q.1: Quality check of raw reads (`00`)

| | |
|---|---|
| Template | `template_00.fastqc.sh` |
| Tool | FastQC 0.12.1 (`qcenv`), MultiQC 1.35 (`multiqcenv`) |
| Input | `input/<dataset>/<s>_1.fq.gz`, `<s>_2.fq.gz` (also `_L1`/`_L2`) |
| Output | `00_fastqc/multiqc_report.html` (all samples on one page, start here), `00_fastqc/multiqc_report_data/` (tables behind it), and per file `00_fastqc/<s>_{1,2}_fastqc.html` and `.zip` (one per read of the pair, so 2 per sample) |
| Config | `step00.dir1` input, `dir2` output, `par1` threads (one file per thread) |
| Notes | Read-only check, nothing downstream depends on it. The report summarises the FastQC modules (pass/warn/fail) for each file. |

## P.1: Adapter and quality trimming (`00-01`)

| | |
|---|---|
| Template | `template_00-01.extract_and_adapter.sh` |
| Tool | AdapterRemoval 2.3.3 (`bioenv`) |
| Input | `input/<dataset>/<s>_1.fq.gz`, `<s>_2.fq.gz` |
| Output | `01_AdapterRemoval/<s>_cleaned_L1.fq.gz`, `<s>_cleaned_L2.fq.gz` (pairs kept), `<s>.singleton.truncated.gz`, `<s>_discarded.fq.gz`, `<s>.settings` (statistics) |
| Config | `step00_01.dir1` input, `dir3` output, `par1` threads, `par2` `--minlength` (50), `par3` `--minquality` (30), `par4` `--maxns` (10), `par5` `--trim5p` (2), `par6` `--trim3p` (2). `dir2`, `dir4`, `dir5` are unused leftovers of the `.tar` input mode. |
| Notes | Adapters are AdapterRemoval's defaults (Illumina TruSeq), `--trimns --trimqualities` enabled. Only paired reads continue; singletons are kept on disk but not used. Output is written with `--gzip-level 1` (fast compression; files ~10-20% larger than the default level, same content). |

## P.2: Host DNA removal (`03`)

| | |
|---|---|
| Template | `template_03.bowtie-remove_host.sh` |
| Tool | Bowtie2 2.3.3.1 (`bowtieenv`), default end-to-end mode |
| Input | `01_AdapterRemoval/<s>_cleaned_L{1,2}.fq.gz`, host index `step03.dir3` |
| Output | `03_bowtie2_output/<s>.unmapped.fq.1.gz`, `<s>.unmapped.fq.2.gz` (pairs that do not align concordantly to the host), `<s>.bowtie2.log` (alignment summary) |
| Config | `step03.dir1` input, `dir2` output, `dir3` index prefix, `par1` threads |
| Notes | The sample name used by every later step becomes `<s>.unmapped`. Alignments are discarded (`-S /dev/null`). |

## C.1 + C.2: Taxonomic profiling and diversity (`04`)

| | |
|---|---|
| Template | `template_04.metaphlan.sh` + `data/calculate_diversity.R` |
| Tool | MetaPhlAn 4.1.0 with database `mpa_vJun23_CHOCOPhlAnSGB_202307` (`/main/db/metaphlan`), R 4.3.3 |
| Input | `03_bowtie2_output/<s>.unmapped.fq.{1,2}.gz` (both mates, comma-separated) |
| Output | `04_metaphlan_output/<s>.unmapped.txt` (profile), `.bowtie2.bz2`, `.sam.bz2`; `merged_abundance_table.txt` (all samples, `merge_metaphlan_tables.py`); `diversity/merged_abundance_table_<metric>.tsv` |
| Config | `step04.dir2` input, `dir3` output, `rscript`, `par1` threads |
| Notes | Alpha metrics: richness, shannon, simpson, gini. Beta: bray-curtis, jaccard, clr (Aitchison). Weighted/unweighted UniFrac are skipped because no `tree.nwk` is provided. The R script installs missing packages at run time (all are already in the image). |

## C.3: Functional profiling (`04b-1`)

| | |
|---|---|
| Template | `template_04b-1.humann_job_array.sh` |
| Tool | HUMAnN 3.9 (`humannenv`), ChocoPhlAn v201901_v31, UniRef90 (DIAMOND) |
| Input | `03_bowtie2_output/<s>.unmapped.fq.1.gz` (**read 1 only**) and the MetaPhlAn profile `04_metaphlan_output/<s>.unmapped.txt` (`--taxonomic-profile`) |
| Output | `04b_humann_output/<s>.unmapped.fq_genefamilies.tsv`, `_pathabundance.tsv`, `_pathcoverage.tsv`, one set per sample (the temporary folder with HUMAnN's own log is deleted, `--remove-temp-output`) |
| Config | `step04b_1.dir1` profiles, `dir2` reads, `dir3` output, `par1` threads, `par2` DIAMOND block size |
| Notes | HUMAnN is pointed to the MetaPhlAn 4.1.0 of `bioenv` (the one that made the profiles): the MetaPhlAn 4.2 pulled into `humannenv` by conda breaks HUMAnN's version check. The profiles were made with database vJun23, which HUMAnN 3.9 maps to ChocoPhlAn via `humann3_sgb_map.Jun23.tsv`. DIAMOND search dominates the run time. `--memory-use maximum` keeps HUMAnN's alignment tables in memory (faster post-processing, same results); `step04b_1.par2` sets DIAMOND's `--block-size` (default 2; larger blocks are faster and use roughly 6 GB of RAM per unit, without changing the hits). Samples whose three tables exist are skipped, so an interrupted step can be resumed. `template_04b-2` is an old SLURM version, not used. |

## C.4: Community-level aggregation (`04b-last`)

| | |
|---|---|
| Template | `template_04b-last.humann_join_tables.sh` |
| Tool | HUMAnN 3.9 utilities (`humann_join_tables`, `humann_renorm_table`, `humann_regroup_table`, `humann_rename_table`, `humann_split_stratified_table`), full utility mapping files |
| Input | the per-sample tables of C.3 |
| Output | `04b_humann_output/merged_tables/`: `humann_joined_{genefamilies,pathabundance,pathcoverage}.tsv`, `_cpm.tsv` versions (copies per million), `merged_table_metacyc_rxn.tsv`, `merged_table_ko.tsv`, `merged_table_ec.tsv`, and `unstratified/` copies without per-species rows |
| Config | `step04b_last.dir1` input, `dir2` output |
| Notes | The stratified tables keep, for each function, the species that carry it; the report uses them to show which species contribute the pathways of each sample. KEGG pathways are not produced (HUMAnN has no KO → pathway mapping); MetaCyc pathways are. The merged taxonomic table is made in C.2 (`04_metaphlan_output/merged_abundance_table.txt`). |

## A.1: Metagenome assembly (`05`)

| | |
|---|---|
| Template | `template_05.spades.assembly.sh` |
| Tool | metaSPAdes 3.13.0 (`spades.py --meta --only-assembler`) |
| Input | `03_bowtie2_output/<s>.unmapped.fq.{1,2}.gz` |
| Output | `05_spades_output/<s>.unmapped/` (SPAdes folder: `contigs.fasta`, `scaffolds.fasta`, `spades.log`, ...) |
| Config | `step05.var1` input, `var2` output, `par1` threads, `par2` memory limit in GB |
| Notes | One assembly per sample (no co-assembly). Error correction is skipped (`--only-assembler`). Memory grows with depth and diversity; SPAdes stops when it reaches `par2`. |

## A.2: Contig filtering (`06`)

| | |
|---|---|
| Template | `template_06.contig_filter.sh` |
| Input | `05_spades_output/*/contigs.fasta` |
| Output | `05_spades_output/contigs/<s>.unmapped.fasta` (copies), `05_spades_output/contigs/filtered/<s>.unmapped.fasta_sort.fasta` (contigs > 500 bp) |
| Config | `step06.var1` SPAdes folder, `var2`, `var3` sub-folder names, `filter1` (unused) |
| Notes | The 500 bp cut-off is hard-coded in the awk command. |

## A.3: Contig index and read mapping (`07`, `08-09`)

| | |
|---|---|
| Templates | `template_07.bowtie-MAG_database.sh`, `template_08-09.mapp_coverage.sh` |
| Tool | Bowtie2 (`--very-sensitive-local`), samtools 1.6 |
| Input | filtered contigs; `03_bowtie2_output/<s>.unmapped.fq.{1,2}.gz` |
| Output | `07_Bowtie_Index/<s>.unmapped.fasta_sort_index_base.*.bt2` (built with `--threads`; note that `-p` in `bowtie2-build` means *packed*, not threads); `08_mapping_coverage/<s>.unmapped.bam`, `.sorted.bam`, `.sorted.bam.bai`, `.sorted.bam.idxstat`, `_bowtie2.log` |
| Config | `step07.dir1` index output, `dir2` contigs, `par1` threads; `step08_09.dir1` reads, `dir2` indexes, `dir3` output, `par1` threads |
| Notes | Each sample's reads are mapped to its own contigs only (single-sample coverage). The unsorted BAM is kept and can be deleted. |

## A.4: Contig coverage and binning (`10-11`)

| | |
|---|---|
| Template | `template_10-11-11b.metabat.sh` |
| Tool | MetaBAT2 2.17 (`jgi_summarize_bam_contig_depths`, `metabat2`) |
| Input | `08_mapping_coverage/*.sorted.bam`, filtered contigs |
| Output | `10_metabat_depth/<s>.unmapped.depth.txt`; `11_metabat_MAG/<s>.unmapped.<n>.fa` (one FASTA per bin) |
| Config | `step10_11.dir1`-`dir5` folders, `par1` minimum contig length (`-m`, 1500), `par2` threads |
| Notes | Mode is hard-coded to 1 (default MetaBAT2 settings into `11_metabat_MAG`); mode 2 (`--minContig 2000 --maxEdges 200` into `11b_metabat_MAG`) is available in the template. |

## A.5: MAG quality (`12-13`)

| | |
|---|---|
| Template | `template_12-13.checkm.sh` |
| Tool | CheckM2 1.0.2 (`checkm2env`, database `uniref100.KO.1.dmnd`); CheckM 1 is installed but not run |
| Input | `11_metabat_MAG/*.fa` |
| Output | `13_checkm2/quality_report.tsv` (Name, Completeness, Contamination, Completeness_Model_Used, Contig_N50, Genome_Size, GC_Content, ...) |
| Config | `step12_13.var1` MAGs, `var2` CheckM output, `var3` CheckM2 output, `par1`/`par2` CheckM threads, `par3` CheckM2 threads |
| Notes | The tool choice is hard-coded to `2` (CheckM2). |

## A.6: MAG filtering and metadata (`14`)

| | |
|---|---|
| Template | `template_14.filter-metadata.sh` |
| Input | `13_checkm2/quality_report.tsv`, `11_metabat_MAG/*.fa` |
| Output | `14_MAGs_high_quality/<MAG>.fa` (MAGs with completeness > `var5` and contamination < `var6`), `my-metadata2.txt` (TORMES metadata) |
| Config | `step14.var1` CheckM2 report, `var2` MAG folder, `var3` output, `var4` folder for mode 2, `var5` completeness (50), `var6` contamination (10), `var7`/`var8` metadata file names |
| Notes | Mode from argument 1, else `STEP14_CHOICE`, else asked: `1` filter (default in the examples), `2` metadata for all MAGs in `var4`, `3` both. The selected MAGs are **moved** out of `11_metabat_MAG/`. If no MAG passes, the step fails and TORMES is not run. 50/10 is the "medium quality" MIMAG threshold. |

## A.7: MAG annotation (`15`)

| | |
|---|---|
| Template | `template_15.tormes_MAG.sh` |
| Tool | TORMES 1.3.0 (`tormes-1.3.0`) |
| Input | `14_MAGs_high_quality/my-metadata2.txt` |
| Output | `15_tormes_MAGs/` (TORMES HTML report, annotation, MLST, resistance and virulence genes, ...) |
| Config | `step15.var1` working dir, `var2` metadata, `par1` output folder, `par2` threads |

## R.1: Integrated report, both branches (`16`)

| | |
|---|---|
| Template | `template_16.report.sh`, which runs `scripts/report/make_report.py` (Python 3 standard library only) |
| Input | the outputs of every step of both branches (dashed purple lines in the figure): Q.1 `00_fastqc/*_fastqc.zip`, P.1.2 `01_AdapterRemoval/*.settings`, P.2.2 `03_bowtie2_output/*.bowtie2.log`, C.1.1 `04_metaphlan_output/<s>.txt`, C.2.1-C.2.2 `merged_abundance_table.txt` and `diversity/`, C.3-C.4 `04b_humann_output/`, A.1-A.2 `05_spades_output/contigs/`, A.8 `06b_quast/transposed_report.tsv`, A.3.2 `08_mapping_coverage/*_bowtie2.log`, A.4.1 `11_metabat_MAG/`, A.4.2 `10_metabat_depth/*.depth.txt`, A.5 `13_checkm2/`, A.6 `14_MAGs_high_quality/`, A.7 `15_tormes_MAGs/` (Kraken2 taxonomy, links to the per-MAG reports); plus `config.yml`, `step_times.tsv` and the run log |
| Output | `report/report.html`: one self-contained interactive file (no internet needed), suitable to send by e-mail |
| Config | `step16.dir1` run folder |
| Notes | One section per step, then a joint section that combines the branches: a per-sample summary (microbial reads, species, Shannon, pathways, contig N50, mapping rate, bins, kept MAGs); the share of each sample's microbial reads that ends up in kept MAGs, other bins, unbinned contigs or is not assembled (contig length × depth from `10_metabat_depth/`, scaled by the A.3 mapping rate); and each MAG's share of the reads next to the MetaPhlAn abundance of the same species and genus in the same sample (MAG taxonomy from TORMES' Kraken2, `15_tormes_MAGs/**/taxonomic-identification-kraken2.txt`; Kraken2 uses the *Lactobacillus* names from before the 2020 genus split, so a MAG is also matched to a MetaPhlAn species of the new genera with the same epithet, e.g. *Lactobacillus rhamnosus* = *Lacticaseibacillus rhamnosus*; if no species matches, the genus value sums *Lactobacillus* and the genera split from it). It also runs when an earlier step fails, so a partial run still gets a report. Regenerate it with `bash scripts/run_pipeline.sh --name <run> --from 16`, or on the host with `python3 scripts/report/make_report.py output/<run>`. |

## A.8: Assembly quality (`06b`)

| | |
|---|---|
| Template | `template_06b.quast.sh` |
| Tool | QUAST 5.3.0 (`quastenv`), no reference genome |
| Input | `05_spades_output/contigs/filtered/*.fasta_sort.fasta` (A.2) |
| Output | `06b_quast/report.html` (all samples side by side), `report.tsv`, `transposed_report.tsv` (one row per sample, read by the report), `icarus.html` (contig size viewer), `basic_stats/` (plots) |
| Config | `step06b.dir1` filtered contigs, `dir2` output folder, `par1` threads, `par2` minimum contig length counted (500) |
| Notes | One QUAST run on all samples, labelled with the sample names. Without a reference QUAST reports contiguity only (contigs, total length, largest contig, N50/N90, L50/L90, GC, Ns), not misassemblies or genome fraction; metaQUAST is not used because it downloads references at run time. Nothing downstream depends on it, so `--from 07` skips it safely. |

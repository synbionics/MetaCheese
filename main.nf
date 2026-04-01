#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================
// 1) PARAMETRI GLOBALI
// =====================
params.input_sample = params.input_sample ?: error("Manca --input_sample (cartella dentro input/)")
params.codice       = params.codice      ?: error("Manca --codice (nome identificativo es. test1)")

// Genera data di oggi (YYYYMMDD)
def today       = new Date().format('yyyyMMdd')
// Crea nome cartella output unico
def full_codice = "${today}_${params.codice}"

// Permetti override da CLI di output_base, altrimenti genera dinamico
params.output_base = "output/${full_codice}"

// =====================
// 2) INCLUDO I MODULI
// =====================

include { extract_tar              } from './modules/00.extract_tar.nf'
include { adapter_removal          } from './modules/01.adapter_removal.nf'
include { bowtie2_map              } from './modules/03.bowtie2_map.nf'

include { metaphlan_sample         } from './modules/04.metaphlan.nf'
include { metaphlan_merge_diversity} from './modules/04.metaphlan_merge_diversity.nf'
include { humann_merge             } from './modules/04b.humann_merge.nf'
include { humann_sample            } from './modules/04b.humann_sample.nf'

include { spades_assembly          } from './modules/05.spades.assembly.nf'
include { contig_filter            } from './modules/06.contig_filter.nf'
include { bowtie_mag_database      } from './modules/07.bowtie_MAG_database.nf'
include { mapp_coverage            } from './modules/08-09.mapp_coverage.nf'
include { metabat_depth            } from './modules/10.metabat_depth.nf'
include { metabat_binning          } from './modules/11-11b.metabat_binning.nf'
include { checkm_lineage_wf        } from './modules/12.checkm_lineage_wf.nf'
include { checkm2_predict          } from './modules/13.checkm2.nf'
include { combine_quality_reports  } from './modules/combine_quality_reports.nf'
include { filter_high_quality_mags } from './modules/14.filter_high_quality_mags.nf'
include { tormes_MAG               } from './modules/15.tormes_MAG.nf'

// =====================
// 3) WORKFLOW PRINCIPALE
// =====================

workflow {

    def config = file(params.config_yaml).text
    def myparams = readYaml(file(params.config_yaml))

    def minlenght  = myparams.adapterRemoval.minlenght
    def minquality = myparams.adapterRemoval.minquality
    def maxns      = myparams.adapterRemoval.maxns
    def trim5p     = myparams.adapterRemoval.trim5p
    def trim3p     = myparams.adapterRemoval.trim3p

    // --------------------------------
    // 3.1) Estrazione archivi .tar → .fq.gz
    // --------------------------------
    Channel
        .fromPath("input/${params.input_sample}/*.tar")
        .set { tar_files }

    // Output: Channel di fq.gz estratti
    def fqgz_files = extract_tar(tar_files).fqgz
    // fqgz_files.view { "DEBUG extract_tar → $it" }

    // --------------------------------
    // 3.2) Raggruppo le coppie di fastq (paired-end)
    // --------------------------------
    def paired_reads = Channel
        .fromFilePairs("${params.output_base}/00_Extract/*_{1,2}.fq.gz")
    // paired_reads.view { "DEBUG paired_reads → $it" }

    // --------------------------------
    // 3.3) Adapter removal
    // --------------------------------
    def cleaned_reads = adapter_removal(paired_reads, minlenght, minquality, maxns, trim5p, trim3p).cleaned
    // cleaned_reads.view { "DEBUG adapter_removal → $it" }

    // --------------------------------
    // 3.4) Host-mapping (Bowtie2) e prendo solo gli unmapped
    // --------------------------------
    def unmapped_reads = bowtie2_map(cleaned_reads).unmapped
    // unmapped_reads.view { "DEBUG bowtie2_map → $it" }

    // --------------------------------
    // 3.5b) Metaphlan profiling
    // --------------------------------
    def metaphlan_results = metaphlan_sample(unmapped_reads).metaphlan_out
    def txts_merged = metaphlan_results.map { it[3] }.collect()
    //txts_merged.view { "DEBUG txts_merged: $it" }
    metaphlan_merge_diversity(txts_merged)

    //unmapped_reads.view { "UNMAPPED: $it" }
    //metaphlan_results.view { "METAPHLAN: $it" }

    def humann_inputs = unmapped_reads.join(metaphlan_results, by: 0)
        .map { it -> tuple(it[0], it[1], it[-1]) }

    //humann_inputs.view { "DEBUG HUMANN INPUT: $it" }

    def humann_outputs = humann_sample(humann_inputs)

    def humann_genefamilies = humann_outputs
        .map { dir -> file("${dir}/*_genefamilies.tsv") }
        .collect()

    humann_merge(humann_genefamilies)

    // --------------------------------
    // 3.5) Assemblaggio SPAdes
    // --------------------------------
    def spades_results = spades_assembly(unmapped_reads).assembly
    // spades_results.view { "DEBUG spades_assembly → $it" }

    // --------------------------------
    // 3.6) Filtraggio contigs
    // --------------------------------
    def filtered_contigs = contig_filter(spades_results).filtered
    // filtered_contigs.view { "DEBUG contig_filter → $it" }

    // --------------------------------
    // 3.7) Creazione indice Bowtie2 per i contig (MAG DB)
    //    Output: Channel di liste [sample, idx1.bt2, idx2.bt2, ...]
    // --------------------------------
    def raw_indexes = bowtie_mag_database(filtered_contigs).db
    // raw_indexes.view { "DEBUG raw_indexes → $it" }

    // --------------------------------
    // 3.8) Ricavo i prefissi degli indici Bowtie2 per ogni sample
    //    Output: tuple(sample, prefix_path)
    // --------------------------------
    def bowtie_indexes = raw_indexes
        .map { record ->
            // record: [sample, idx1.bt2, ...]
            def sample   = record[0]
            def bt2_files = record[1..-1]
            def idx1 = bt2_files.find { it.name.endsWith('.1.bt2') }
            if (!idx1) error "Nessun file '*.1.bt2' per il sample ${sample}"
            def prefix = idx1.toString().replaceAll(/\.1\.bt2$/, '')
            tuple(sample, prefix)
        }
    // bowtie_indexes.view { "DEBUG bowtie_indexes → $it" }

    // --------------------------------
    // 3.9) Preparo input per mapp_coverage: [sample, fq1, fq2, idx_prefix]
    // --------------------------------
    def mapping_inputs = unmapped_reads
        .join(bowtie_indexes, by: 0)
        .map { sample, fq1, fq2, idx_prefix -> tuple(sample, fq1, fq2, idx_prefix) }
    // mapping_inputs.view { "DEBUG mapping_inputs → $it" }

    def mapp_coverage_output = mapp_coverage(mapping_inputs)

    // --------------------------------
    // 3.10) Calcolo profondità di copertura (jgi_summarize_bam_contig_depths)
    //      Prende [sample, sorted_bam] in input
    // --------------------------------
    def sorted_bams = mapp_coverage_output.map { t -> tuple(t[0], t[2]) }
    def depth_files = metabat_depth(sorted_bams)

    // --------------------------------
    // 3.11) Binning con MetaBAT2
    // --------------------------------
    def metabat_join = filtered_contigs.join(depth_files)
    def metabat_bins = metabat_binning(metabat_join)

    // --------------------------------
    // 3.12) CheckM lineage_wf sui bins (batch)
    // --------------------------------
    checkm_lineage_wf(metabat_bins)

    // --------------------------------
    // 3.13) CheckM2 predict sui bins (batch)
    //      Prepara input raggruppando i file per sample
    // --------------------------------
    def checkm2_inputs = metabat_bins
        .groupTuple(by: 0)
        .map { sample, files -> tuple(sample, files.flatten()) }

    //checkm2_inputs.view { "DEBUG checkm2_inputs → $it" }


    def checkm2_outputs = checkm2_predict(checkm2_inputs)

    // raccogliamo solo i path dei .tsv
    def all_checkm2_outputs = checkm2_outputs.map { it[1] }.collect()

    def combined_report = combine_quality_reports(all_checkm2_outputs)
    //combined_report.view { "DEBUG combined_report → $it" }

    // --------------------------------
    // 3.14 - 15) Filtro i MAGs di alta qualità
    // ---------------------------------

    def mag_files = metabat_bins.map { it[1] }.flatten()
    def all_mags = mag_files.collect()

    def filtered_mags = filter_high_quality_mags(combined_report, all_mags)

}

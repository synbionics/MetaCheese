import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "regression"))
import summarize as su  # noqa: E402


def write(run, rel, text):
    path = os.path.join(run, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def tiny_run(run):
    write(run, "01_AdapterRemoval/S1.settings", "Total number of read pairs: 100\nNumber of retained reads: 190\n")
    write(run, "03_bowtie2_output/S1.bowtie2.log",
          "100 reads; of these:\n  90 (90.00%) aligned concordantly 0 times\n11.00% overall alignment rate\n")
    write(run, "04_metaphlan_output/merged_abundance_table.txt",
          "#mpa_db\nclade_name\tS1.unmapped\nk__B|p__P|c__C|o__O|f__F|g__G|s__G_sp\t100.0\n"
          "k__B|p__P|c__C|o__O|f__F|g__G|s__G_absent\t0.0\n")
    write(run, "04_metaphlan_output/diversity/merged_abundance_table_shannon.tsv", "diversity_shannon\nS1.unmapped\t1.5\n")
    write(run, "04b_humann_output/S1.unmapped.fq_pathabundance.tsv",
          "# Pathway\tS1\nUNMAPPED\t25\nUNINTEGRATED\t25\nPWY-1\t50\nPWY-1|g__G.s__G_sp\t50\nPWY-0\t0\n")
    write(run, "04b_humann_output/S1.unmapped.fq_genefamilies.tsv",
          "# Gene\tS1\nUNMAPPED\t5\nU1\t1\nU1|g__G.s__G_sp\t1\nU2\t1\n")
    write(run, "13_checkm2/quality_report.tsv",
          "Name\tCompleteness\tContamination\tGenome_Size\nS1.unmapped.1\t90\t2\t1000\nS1.unmapped.2\t60\t4\t1000\nS1.unmapped.3\t20\t1\t1000\n")
    write(run, "14_MAGs_high_quality/S1.unmapped.1.fa", ">c\nACGT\n")
    write(run, "14_MAGs_high_quality/S1.unmapped.2.fa", ">d\nACGT\n")
    write(run, "15_tormes_MAGs/taxonomic_identification/taxonomic-identification-kraken2.txt",
          "Isolate\tFamily\tGenus\tSpecies\nS1.unmapped.1\tF (100%)\tG (100%)\tG sp (99.0%)\n")
    write(run, "step_times.tsv", "00.fastqc.sh\t3\n17.report.sh\t4\n")
    write(run, "versions.tsv", "Tool\t1.0\n")


class Summarize(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = self.tmp.name
        tiny_run(self.run)
        self.summary = su.summarize(self.run, 0)
        self.s1 = self.summary["per_sample"]["S1"]

    def tearDown(self):
        self.tmp.cleanup()

    def test_reads(self):
        self.assertEqual(self.summary["samples"], ["S1"])
        self.assertEqual((self.s1["trim_pairs"], self.s1["microbial_pairs"], self.s1["host_pairs"]), (100, 90, 10))

    def test_only_species_found_in_the_sample_are_kept(self):
        self.assertEqual(self.s1["species"], {"G sp": 100.0})

    def test_humann_ignores_stratified_lines(self):
        self.assertEqual(self.s1["pathways"], {"PWY-1": 50.0})
        self.assertEqual(self.s1["humann_unmapped_pct"], 25.0)
        self.assertEqual(self.s1["humann_unintegrated_pct"], 25.0)
        self.assertEqual(self.s1["gene_families"], 2)

    def test_mags(self):
        self.assertEqual((self.s1["mags_total"], self.s1["mags_kept"]), (3, 2))
        self.assertEqual(self.s1["mag_mean_completeness"], 75.0)
        self.assertEqual(self.s1["mag_species"], ["G sp", "unclassified"])

    def test_joint_values_need_step_16(self):
        self.assertIsNone(self.s1["kept_read_share"])
        self.assertIsNone(self.s1["kept_mags_metaphlan_species"])
        write(self.run, "16_integration/read_fate.tsv",
              "sample\tkept_mags_pct\tother_bins_pct\tunbinned_pct\tnot_assembled_pct\nS1\t40\t10\t30\t20\n")
        write(self.run, "16_integration/mags.tsv",
              "mag\tsample\tkept\tcompleteness\tcontamination\tread_share_pct\tkraken2_species\tkraken2_genus\t"
              "metaphlan_name\tmetaphlan_species_pct\tmetaphlan_genus_pct\n"
              "S1.unmapped.1\tS1\tyes\t90\t2\t30\tG sp\tG\tG sp\t100\t100\n"
              "S1.unmapped.2\tS1\tyes\t60\t4\t10\t\t\t\t\t\n")
        s1 = su.summarize(self.run, 0)["per_sample"]["S1"]
        self.assertEqual((s1["kept_read_share"], s1["kept_mags_metaphlan_species"]), (40.0, 1))

    def test_files_and_steps(self):
        files = self.summary["files"]
        self.assertTrue(files["01_AdapterRemoval/S1.settings"])
        self.assertFalse(files["08_mapping_coverage/S1.unmapped.sorted.bam"])
        self.assertEqual(self.summary["steps"], ["00.fastqc.sh", "17.report.sh"])
        self.assertEqual(self.summary["versions"], {"Tool": "1.0"})


if __name__ == "__main__":
    unittest.main()

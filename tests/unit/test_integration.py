import os
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts", "integration"))
sys.path.insert(0, os.path.join(REPO, "scripts", "report"))
import integrate as ig  # noqa: E402
import make_report as mr  # noqa: E402


def write(run, rel, text):
    path = os.path.join(run, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def mag(name, sample="S1", passed=True):
    return {"name": name, "sample": sample, "passed": passed, "completeness": 90.0, "contamination": 1.0}


def profile(*rows):
    """MetaPhlAn profile of one sample S1: (rank, taxon, %) rows."""
    out = {"samples": ["S1"], "profiles": {}}
    for rank, taxon, value in rows:
        out["profiles"].setdefault(rank, []).append({"taxon": taxon, "values": [value]})
    return out


class Integration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def kraken(self, *rows):
        write(self.run, "15_tormes_MAGs/taxonomic_identification/taxonomic-identification-kraken2.txt",
              "Isolate\tFamily\tGenus\tSpecies\n" + "".join(f"{m}\tF (100%)\t{g} (100%)\t{s} (99%)\n" for m, g, s in rows))

    def test_read_fate_splits_mapped_reads_by_contig_fate(self):
        write(self.run, "11_metabat_MAG/S1.unmapped.1.fa", ">c1\nA\n")
        write(self.run, "11_metabat_MAG/S1.unmapped.2.fa", ">c2\nA\n")
        write(self.run, "14_MAGs_high_quality/S1.unmapped.1.fa", ">c1\nA\n")
        # length x depth: c1 400 (kept MAG), c2 200 (other bin), c3 400 (unbinned)
        write(self.run, "10_metabat_depth/S1.unmapped.depth.txt",
              "contigName\tcontigLen\ttotalAvgDepth\nc1\t100\t4\nc2\t100\t2\nc3\t200\t2\n")
        fate, share = ig.read_fate(self.run, ["S1"], {"S1": 80.0})
        self.assertEqual(fate["S1"], {"kept": 32.0, "other_bins": 16.0, "unbinned": 32.0, "not_assembled": 20.0})
        self.assertEqual(share, {"S1.unmapped.1": 32.0, "S1.unmapped.2": 16.0})

    def test_species_found_in_the_same_sample(self):
        self.kraken(("S1.unmapped.1", "Streptococcus", "Streptococcus thermophilus"))
        prof = profile(("species", "Streptococcus thermophilus", 40.0), ("genus", "Streptococcus", 41.0))
        (m,), notes = ig.mags_vs_metaphlan(self.run, [mag("S1.unmapped.1")], prof, {"S1.unmapped.1": 30.0})
        self.assertEqual((m["mpa_name"], m["mpa_species"], m["mpa_genus"], m["share"]),
                         ("Streptococcus thermophilus", 40.0, 41.0, 30.0))
        self.assertEqual(notes, [])

    def test_absent_species_is_not_matched(self):
        self.kraken(("S1.unmapped.1", "Streptococcus", "Streptococcus thermophilus"))
        prof = profile(("species", "Streptococcus thermophilus", 0.0), ("genus", "Streptococcus", 2.0))
        (m,), _ = ig.mags_vs_metaphlan(self.run, [mag("S1.unmapped.1")], prof, {})
        self.assertEqual((m["mpa_name"], m["mpa_species"], m["mpa_genus"]), (None, None, 2.0))

    def test_old_lactobacillus_name_matches_new_genus(self):
        self.kraken(("S1.unmapped.1", "Lactobacillus", "Lactobacillus rhamnosus"))
        prof = profile(("species", "Lacticaseibacillus rhamnosus", 50.0), ("genus", "Lacticaseibacillus", 50.0))
        (m,), notes = ig.mags_vs_metaphlan(self.run, [mag("S1.unmapped.1")], prof, {})
        self.assertEqual((m["mpa_name"], m["mpa_species"], m["mpa_genus"]), ("Lacticaseibacillus rhamnosus", 50.0, 50.0))
        self.assertEqual(notes, [ig.RENAMED_GENERA["lactobacillus"]["note"]])

    def test_unmatched_lactobacillus_sums_the_split_genera(self):
        self.kraken(("S1.unmapped.1", "Lactobacillus", "Lactobacillus sp"))
        prof = profile(("genus", "Lacticaseibacillus", 30.0), ("genus", "Lactobacillus", 5.0),
                       ("genus", "Streptococcus", 60.0))
        (m,), notes = ig.mags_vs_metaphlan(self.run, [mag("S1.unmapped.1")], prof, {})
        self.assertEqual((m["mpa_name"], m["mpa_genus"]), (None, 35.0))
        self.assertEqual(notes, [ig.RENAMED_GENERA["lactobacillus"]["genus_note"]])

    def test_report_reads_back_the_tables(self):
        self.kraken(("S1.unmapped.1", "Streptococcus", "Streptococcus thermophilus"))
        write(self.run, "11_metabat_MAG/S1.unmapped.1.fa", ">c1\nA\n")
        write(self.run, "14_MAGs_high_quality/S1.unmapped.1.fa", ">c1\nA\n")
        write(self.run, "10_metabat_depth/S1.unmapped.depth.txt", "contigName\tcontigLen\ttotalAvgDepth\nc1\t100\t4\n")
        write(self.run, "08_mapping_coverage/S1.unmapped_bowtie2.log",
              "10 reads; of these:\n  2 (20.00%) aligned concordantly 0 times\n90.00% overall alignment rate\n")
        write(self.run, "04_metaphlan_output/merged_abundance_table.txt",
              "#mpa\nclade_name\tS1.unmapped\nk__B|p__P|c__C|o__O|f__F|g__Streptococcus\t70\n"
              "k__B|p__P|c__C|o__O|f__F|g__Streptococcus|s__Streptococcus_thermophilus\t70\n")
        write(self.run, "01_AdapterRemoval/S1.settings", "Total number of read pairs: 10\n")
        argv, sys.argv = sys.argv, ["integrate.py", self.run]
        try:
            ig.main()
        finally:
            sys.argv = argv
        joint = mr.collect_joint(self.run)
        self.assertEqual(joint["fate"]["S1"], {"kept": 90.0, "other_bins": 0.0, "unbinned": 0.0, "not_assembled": 10.0})
        (m,) = joint["mags"]
        self.assertEqual((m["name"], m["passed"], m["share"], m["mpa_name"], m["mpa_species"]),
                         ("S1.unmapped.1", True, 90.0, "Streptococcus thermophilus", 70.0))
        self.assertTrue(joint["taxonomy"])
        self.assertEqual(joint["name_notes"], [])

    def functions_run(self):
        """One kept MAG of Streptococcus thermophilus in S1, with HUMAnN tables."""
        write(self.run, "15_tormes_MAGs/annotation/S1.unmapped.1_annotation/S1.unmapped.1.tsv",
              "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct\n"
              "a\tCDS\t9\t\t1.1.1.1\t\tx\nb\tCDS\t9\t\t2.2.2.2\t\tx\nc\tCDS\t9\t\t3.6.3.-\t\tx\n"
              "d\tCDS\t9\t\t\t\tx\ne\tCDS\t9\t\t1.1.1.1\t\tx\n")
        write(self.run, "04b_humann_output/S1.unmapped.fq_pathabundance.tsv",
              "# Pathway\tS1\nUNMAPPED\t10\nUNINTEGRATED\t50\nUNINTEGRATED|g__Streptococcus.s__Streptococcus_thermophilus\t30\n"
              "PWY-1: one\t40\nPWY-1: one|g__Streptococcus.s__Streptococcus_thermophilus\t30\nPWY-1: one|unclassified\t10\n"
              "PWY-2: two\t20\nPWY-2: two|g__Other.s__Other_sp\t20\n")
        write(self.run, "04b_humann_output/S1.unmapped.fq_genefamilies.tsv",
              "# Gene Family\tS1\nUNMAPPED\t50\nU1\t30\nU1|g__Streptococcus.s__Streptococcus_thermophilus\t30\n"
              "U2\t20\nU2|g__Other.s__Other_sp\t20\n")
        write(self.run, "04b_humann_output/merged_tables/merged_table_ec.tsv",
              "# Gene Family\tS1.unmapped.fq_Abundance-CPM\nUNMAPPED\t100\nUNGROUPED\t500\n"
              "1.1.1.1: A\t300\n1.1.1.1: A|g__Streptococcus.s__Streptococcus_thermophilus\t300\n"
              "3.3.3.3: C\t50\n3.3.3.3: C|g__Streptococcus.s__Streptococcus_thermophilus\t50\n"
              "4.4.4.4: D\t50\n4.4.4.4: D|g__Other.s__Other_sp\t50\n")
        joint = [dict(mag("S1.unmapped.1"), species="Streptococcus thermophilus"),
                 dict(mag("S1.unmapped.2"), species="Lactobacillus kefiri"),
                 dict(mag("S1.unmapped.3", passed=False), species=None)]
        return ig.mag_functions(self.run, ["S1"], joint)

    def test_mag_linked_to_humann_species(self):
        summary, pathways, _, _ = self.functions_run()
        st = summary[0]
        self.assertEqual(len(summary), 2)  # kept MAGs only
        self.assertEqual(st["humann_species"], "Streptococcus thermophilus")
        self.assertEqual((st["humann_pathways"], st["humann_pathway_share_pct"]), (1, 50.0))  # 30 of 60 stratified
        self.assertEqual((st["humann_gene_families"], st["humann_gene_family_share_pct"]), (1, 30.0))  # 30 of 100
        self.assertEqual(pathways, [["S1.unmapped.1", "S1", "PWY-1", "one", 30.0, 75.0]])

    def test_mag_without_humann_species(self):
        summary, _, ecs, _ = self.functions_run()
        k = summary[1]
        self.assertIsNone(k["humann_species"])
        self.assertNotIn("humann_pathways", k)
        self.assertNotIn("S1.unmapped.2", [r[0] for r in ecs])  # not annotated, not linked

    def test_ec_comparison(self):
        summary, _, ecs, coverage = self.functions_run()
        st = summary[0]
        # MAG: 1.1.1.1, 2.2.2.2 (partial 3.6.3.- ignored); HUMAnN species: 1.1.1.1, 3.3.3.3
        self.assertEqual((st["mag_ecs"], st["humann_species_ecs"], st["ecs_in_both"],
                          st["ecs_humann_species_only"], st["ecs_mag_only"], st["mag_ecs_in_community"]), (2, 2, 1, 1, 1, 1))
        self.assertEqual(ecs, [["S1.unmapped.1", "S1", "1.1.1.1", "A", True, 300.0, 300.0],
                               ["S1.unmapped.1", "S1", "2.2.2.2", "", True, 0.0, 0.0],
                               ["S1.unmapped.1", "S1", "3.3.3.3", "C", False, 50.0, 50.0]])
        # community ECs 1.1.1.1 (300), 3.3.3.3 (50), 4.4.4.4 (50): the MAG encodes 300 of 400
        self.assertEqual(coverage, [["S1", 3, 1, 75.0]])

    def test_humann_link_across_the_lactobacillus_split(self):
        old = {"lactobacillus rhamnosus", "streptococcus thermophilus"}
        self.assertEqual(ig.humann_key("lacticaseibacillus rhamnosus", old), "lactobacillus rhamnosus")
        self.assertEqual(ig.humann_key("streptococcus thermophilus", old), "streptococcus thermophilus")
        self.assertIsNone(ig.humann_key("streptococcus rhamnosus", old))  # outside the split: names must agree
        new = {"lacticaseibacillus rhamnosus"}
        self.assertEqual(ig.humann_key("lactobacillus rhamnosus", new), "lacticaseibacillus rhamnosus")
        # two candidates with the same epithet: no guess
        self.assertIsNone(ig.humann_key("lactobacillus casei", {"lacticaseibacillus casei", "lactobacillus casei group"} | {"loigolactobacillus casei"}))

    def test_report_without_step_16(self):
        self.assertIsNone(mr.collect_joint(self.run))


if __name__ == "__main__":
    unittest.main()

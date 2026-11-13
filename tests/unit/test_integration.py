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

    def test_report_without_step_16(self):
        self.assertIsNone(mr.collect_joint(self.run))


if __name__ == "__main__":
    unittest.main()

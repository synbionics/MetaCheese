import os
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts", "report"))
import make_report as mr  # noqa: E402

BOWTIE_LOG = """1000 reads; of these:
  1000 (100.00%) were paired; of these:
    900 (90.00%) aligned concordantly 0 times
    90 (9.00%) aligned concordantly exactly 1 time
    10 (1.00%) aligned concordantly >1 times
11.50% overall alignment rate
"""


class Parsers(unittest.TestCase):
    def test_parse_bowtie_log(self):
        self.assertEqual(mr.parse_bt2(BOWTIE_LOG), {"pairs_in": 1000, "not_concordant": 900, "overall_rate": 11.5})
        self.assertIsNone(mr.parse_bt2("not a log"))

    def test_clean_sample(self):
        self.assertEqual(mr.clean_sample("/x/D_GP_02.unmapped.fq_Abundance-CPM"), "D_GP_02")

    def test_match_sample_keeps_similar_names_apart(self):
        samples = ["D_GP_10", "ND_GP_10", "ND_GP_1"]
        self.assertEqual(mr.match_sample("ND_GP_10.unmapped.3", samples), "ND_GP_10")
        self.assertEqual(mr.match_sample("D_GP_10.unmapped", samples), "D_GP_10")
        self.assertEqual(mr.match_sample("ND_GP_1_bowtie2.log", samples), "ND_GP_1")

    def test_n50(self):
        self.assertEqual(mr.n50([]), 0)
        self.assertEqual(mr.n50([100, 50, 25, 25]), 100)
        self.assertEqual(mr.n50([40, 30, 20, 10]), 30)

    def test_host_name_from_any_path(self):
        cfg = 'step03:\n  dir1: "a"\n  dir3: "{}"   # comment\n  par1: 1\n'
        self.assertEqual(mr.host_name(cfg.format("../../data/gene/Bos_taurus/Bos_taurus")), "Bos taurus")
        self.assertEqual(mr.host_name(cfg.format("/data/genomes/Homo_sapiens/GRCh38")), "Homo sapiens")

    def test_collect_trimming(self):
        with tempfile.TemporaryDirectory() as run:
            os.mkdir(os.path.join(run, "01_AdapterRemoval"))
            with open(os.path.join(run, "01_AdapterRemoval", "S1.settings"), "w") as fh:
                fh.write("Total number of read pairs: 100\nNumber of retained reads: 190\n"
                         "Number of reads with adapters[1]: 7\nAverage length of retained reads: 148.5\n")
            got = mr.collect_trimming(run)["S1"]
        self.assertEqual((got["pairs"], got["retained_reads"], got["adapters"], got["avg_length"]), (100, 190, 7, 148.5))

    def test_collect_tormes_puts_the_summary_report_apart(self):
        with tempfile.TemporaryDirectory() as run:
            for rel in ("15_tormes_MAGs/tormes_report.html", "15_tormes_MAGs/report_files/tormes_report.html",
                        "15_tormes_MAGs/genome_stats/M1_genome_stats/report.html"):
                os.makedirs(os.path.dirname(os.path.join(run, rel)), exist_ok=True)
                open(os.path.join(run, rel), "w").close()
            got = mr.collect_tormes(run)
        self.assertEqual(got["summary"], "15_tormes_MAGs/tormes_report.html")
        self.assertEqual(got["reports"], ["15_tormes_MAGs/genome_stats/M1_genome_stats/report.html"])
        self.assertFalse(got["summary_failed"])

    def test_collect_tormes_notes_a_failed_summary_report(self):
        with tempfile.TemporaryDirectory() as run:
            os.mkdir(os.path.join(run, "15_tormes_MAGs"))
            with open(os.path.join(run, "15_tormes_MAGs", "tormes.log"), "w") as fh:
                fh.write("Tormes report started\nWARNING: html report file could not be created\n")
            got = mr.collect_tormes(run)
        self.assertEqual((got["summary"], got["summary_failed"], got["reports"]), (None, True, []))


if __name__ == "__main__":
    unittest.main()

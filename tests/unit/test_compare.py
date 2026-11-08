import copy
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "regression"))
import compare as c  # noqa: E402


def summary():
    sample = {"trim_pairs": 100, "retained_reads": 190, "microbial_pairs": 90, "host_pairs": 5,
              "shannon": 1.5, "richness": 10, "species": {"A a": 60.0, "B b": 40.0},
              "pathways": {"P1": 1.0, "P2": 3.0}, "pathway_count": 2, "gene_families": 50,
              "humann_unmapped_pct": 10.0, "humann_unintegrated_pct": 20.0, "contigs": 500,
              "assembly_bp": 1000000, "n50": 2000, "mapping_rate": 80.0, "mags_total": 3, "mags_kept": 2,
              "mag_mean_completeness": 90.0, "mag_mean_contamination": 1.0, "kept_read_share": 30.0,
              "mag_species": ["A a", "B b"]}
    return {"run": "r", "samples": ["S1"], "steps": ["00", "16"], "files": {"x/S1.txt": True},
            "columns": {"t.tsv": ["a", "b"]}, "versions": {"Tool": "1.0"}, "pipeline_exit": 0,
            "per_sample": {"S1": sample}, "between_samples": {}}


LIMITS = {"microbial_pairs": 0, "shannon": 0.1, "species_bray_curtis": 0.01, "n50": 0.05}


def statuses(new):
    return [(r["status"], r["metric"]) for r in c.compare(summary(), new, LIMITS) if r["status"] != "PASS"]


class Compare(unittest.TestCase):
    def test_identical_runs_pass(self):
        self.assertEqual(statuses(summary()), [])

    def test_failed_pipeline_fails(self):
        new = summary()
        new["pipeline_exit"] = 1
        self.assertIn(("FAIL", "pipeline"), statuses(new))

    def test_missing_sample_fails(self):
        new = summary()
        new["samples"] = []
        self.assertIn(("FAIL", "sample missing"), statuses(new))

    def test_missing_file_fails(self):
        new = summary()
        new["files"]["x/S1.txt"] = False
        self.assertEqual(statuses(new), [("FAIL", "file missing or empty")])

    def test_missing_step_fails(self):
        new = summary()
        new["steps"] = ["00"]
        self.assertEqual(statuses(new), [("FAIL", "step not completed")])

    def test_changed_columns_fail(self):
        new = summary()
        new["columns"]["t.tsv"] = ["a"]
        self.assertEqual(statuses(new), [("FAIL", "table columns changed")])

    def test_missing_value_fails(self):
        new = summary()
        new["per_sample"]["S1"]["n50"] = None
        self.assertEqual(statuses(new), [("FAIL", "n50 missing")])

    def test_value_within_limit_passes(self):
        new = summary()
        new["per_sample"]["S1"]["shannon"] = 1.55
        self.assertEqual(statuses(new), [])

    def test_value_beyond_limit_is_checked(self):
        new = summary()
        new["per_sample"]["S1"]["shannon"] = 1.8
        self.assertEqual(statuses(new), [("CHECK", "shannon")])

    def test_relative_limit(self):
        new = summary()
        new["per_sample"]["S1"]["n50"] = 1950  # 2.5% lower
        self.assertEqual(statuses(new), [])
        new["per_sample"]["S1"]["n50"] = 1000  # 50% lower
        self.assertEqual(statuses(new), [("CHECK", "n50")])

    def test_zero_limit_catches_one_read(self):
        new = summary()
        new["per_sample"]["S1"]["microbial_pairs"] = 91
        self.assertEqual(statuses(new), [("CHECK", "microbial_pairs")])

    def test_changed_version_is_checked(self):
        new = summary()
        new["versions"]["Tool"] = "2.0"
        self.assertEqual(statuses(new), [("CHECK", "version changed")])

    def test_major_species_change_is_checked(self):
        new = summary()
        new["per_sample"]["S1"]["species"] = {"A a": 60.0, "C c": 40.0}
        found = statuses(new)
        self.assertIn(("CHECK", "species lost"), found)
        self.assertIn(("CHECK", "species gained"), found)

    def test_rare_species_change_is_not_listed(self):
        base, new = summary(), summary()
        new["per_sample"]["S1"]["species"]["Rare r"] = 0.5
        metrics = [r["metric"] for r in c.compare(base, new, {}) if r["status"] != "PASS"]
        self.assertNotIn("species gained", metrics)

    def test_mag_species_compared_as_a_multiset(self):
        new = summary()
        new["per_sample"]["S1"]["mag_species"] = ["A a", "A a"]
        found = statuses(new)
        self.assertIn(("CHECK", "MAG species lost"), found)
        self.assertIn(("CHECK", "MAG species gained"), found)


class Helpers(unittest.TestCase):
    def test_bray_curtis(self):
        self.assertEqual(c.bray_curtis({"a": 1}, {"a": 5}), 0)
        self.assertEqual(c.bray_curtis({"a": 1}, {"b": 1}), 1)
        self.assertAlmostEqual(c.bray_curtis({"a": 1, "b": 1}, {"a": 1}), 0.5)
        self.assertEqual(c.bray_curtis({}, {}), 0)
        self.assertEqual(c.bray_curtis({}, {"a": 1}), 1)

    def test_suggested_limits_cover_the_noise(self):
        new = summary()
        new["per_sample"]["S1"]["shannon"] = 1.6
        rows = c.compare(summary(), new, {})
        limits = c.suggest_limits(rows)
        self.assertGreaterEqual(limits["shannon"], 0.2)
        self.assertEqual(limits["microbial_pairs"], 0)
        again = [r for r in c.compare(summary(), new, limits) if r["status"] != "PASS"]
        self.assertEqual(again, [])


if __name__ == "__main__":
    unittest.main()

import os
import re
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def read(*parts):
    with open(os.path.join(REPO, *parts), encoding="utf-8") as fh:
        return fh.read()


def host_setting(text):
    return re.search(r'^step03:\s*\n(?:[ \t]+.*\n)*?[ \t]+dir3:\s*"([^"]*)"', text, re.M).group(1)


class Config(unittest.TestCase):
    def test_general_config_has_no_default_host(self):
        self.assertEqual(host_setting(read("config", "Config.yml")), "")

    def test_examples_set_their_host(self):
        for example in ("PRJNA1198075", "small_PRJNA1198075"):
            self.assertEqual(host_setting(read("examples", example, "Config.yml")),
                             "../../data/gene/Bos_taurus/Bos_taurus")

    def test_configurator_embeds_the_current_config(self):
        html = read("config", "configurator.html")
        for line in read("config", "Config.yml").splitlines():
            if line.strip():
                self.assertIn(line, html)

    def test_general_code_names_no_host_or_taxon_default(self):
        for rel in ("examples/run_in_docker.sh", "scripts/run_pipeline.sh", "scripts/build_bowtie2_index.sh"):
            code = [l for l in read(*rel.split("/")).splitlines() if not l.lstrip().startswith("#")]
            self.assertNotIn("Bos_taurus", "\n".join(code), rel)


if __name__ == "__main__":
    unittest.main()

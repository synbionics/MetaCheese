"""run_in_docker.sh host checks, with a stand-in for docker so no image is needed."""
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUN_NAME = "unit_test_run_in_docker"
FAKE_DOCKER = """#!/bin/sh
[ "$1" = image ] && exit 0
printf '%s\\n' "$@" > "$DOCKER_ARGS"
"""


class Wrapper(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        os.mkdir(os.path.join(self.tmp, "bin"))
        docker = os.path.join(self.tmp, "bin", "docker")
        with open(docker, "w") as fh:
            fh.write(FAKE_DOCKER)
        os.chmod(docker, os.stat(docker).st_mode | stat.S_IEXEC)
        self.args_file = os.path.join(self.tmp, "args")
        self.host = os.path.join(self.tmp, "host")
        os.mkdir(self.host)

    def tearDown(self):
        shutil.rmtree(self.tmp)
        shutil.rmtree(os.path.join(REPO, "output", RUN_NAME), ignore_errors=True)

    def run_wrapper(self, dir3):
        config = os.path.join(self.tmp, "config.yml")
        with open(config, "w") as fh:
            fh.write(f'step00:\n  dir3: "other"\nstep03:\n  dir1: "x"\n  dir3: "{dir3}"   # comment\n  par1: 1\nstep04:\n  dir3: "other"\n')
        env = {**os.environ, "PATH": f"{self.tmp}/bin:{os.environ['PATH']}", "DOCKER_ARGS": self.args_file}
        return subprocess.run(["bash", os.path.join(REPO, "examples", "run_in_docker.sh"), "--name", RUN_NAME, "--config", config],
                              env=env, capture_output=True, text=True)

    def docker_args(self):
        with open(self.args_file) as fh:
            return fh.read().splitlines()

    def assert_rejected(self, result, message):
        self.assertEqual(result.returncode, 1)
        self.assertIn(message, result.stderr)
        self.assertFalse(os.path.exists(os.path.join(REPO, "output", RUN_NAME)), "no run folder for a rejected run")

    def test_no_host_is_rejected(self):
        self.assert_rejected(self.run_wrapper(""), "host removal needs a host genome")

    def test_relative_path_outside_data_is_rejected(self):
        self.assert_rejected(self.run_wrapper("genomes/x"), "host removal needs a host genome")

    def test_path_with_spaces_is_rejected(self):
        self.assert_rejected(self.run_wrapper("/data/my genomes/x"), "spaces")

    def test_folder_without_genome_is_rejected(self):
        self.assert_rejected(self.run_wrapper(f"{self.host}/x"), "no host genome")

    def test_fasta_in_any_folder_is_mounted_writable(self):
        open(os.path.join(self.host, "g.fna"), "w").close()
        result = self.run_wrapper(f"{self.host}/x")
        self.assertEqual(result.returncode, 0, result.stderr)
        args = self.docker_args()
        self.assertIn(f"{self.host}:{self.host}", args)
        self.assertIn(f"HOST_INDEX={self.host}/x", args)

    def test_existing_index_is_mounted_read_only(self):
        open(os.path.join(self.host, "x.1.bt2"), "w").close()
        result = self.run_wrapper(f"{self.host}/x")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"{self.host}:{self.host}:ro", self.docker_args())

    def test_saved_config_wins_when_resuming(self):
        run = os.path.join(REPO, "output", RUN_NAME)
        os.makedirs(run)
        open(os.path.join(self.host, "saved.1.bt2"), "w").close()
        with open(os.path.join(run, "config.yml"), "w") as fh:
            fh.write(f'step03:\n  dir3: "{self.host}/saved"\n')
        result = self.run_wrapper("")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"HOST_INDEX={self.host}/saved", self.docker_args())


if __name__ == "__main__":
    unittest.main()

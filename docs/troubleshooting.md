# Troubleshooting

**The run stopped (error, reboot, power cut).** Run the same command again
with the same `--name` and `--from` the step that stopped, e.g.
`bash examples/run_in_docker.sh --input my_dataset --name my_run --config config/Config.yml --from 04b-1`.
The run reuses its saved `output/<run>/config.yml`. HUMAnN skips the samples
that are already finished and removes the half-written temporary folder of the
interrupted one.

**`Permission denied` on files or folders in `output/<run>/`.** They were
written as root by the container and the run ended before
`examples/run_in_docker.sh` could give them back to your user (e.g. after a
reboot). Resuming the run fixes it, or run
`sudo chown -R $USER output/<run>`.

**`ERROR: key '<step>.<key>' missing`.** The run folder was created with an
older config. Copy the key from `config/Config.yml` into
`output/<run>/config.yml`.

**metaSPAdes killed / out of memory.** Raise `step05.par2` on a larger
machine, or run fewer samples at a time. SPAdes memory grows with read depth
and community complexity.

**No MAG passes step 14, TORMES does not run.** Expected with shallow data:
with few reads per sample the bins are small and incomplete. Use deeper data,
or lower `step14.var5` (completeness) for testing only.

**`Got permission denied while trying to connect to the Docker daemon`.**
Add yourself to the docker group, `sudo usermod -aG docker $USER`, and log in
again.

**`yq: command not found`.** `scripts/run_pipeline.sh` must run inside the
container, where `yq` is installed. From the host use
`examples/run_in_docker.sh`.

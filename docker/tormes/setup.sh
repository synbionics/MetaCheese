#!/usr/bin/env bash
# Pinned replacement for TORMES 1.3.0's own `tormes-setup` (run at image build).
#
# tormes-setup clones git repositories and downloads files from `master`
# branches, and installs Bioconductor packages with BiocManager, so two builds
# can differ. This script does the same steps at the versions of the image that
# passed the regression test (2026-11), and leaves out two of them:
# - the minikraken_8GB_202003 download: the Dockerfile installs the Kraken2
#   Standard-8 2026-10-22 database instead;
# - the Kaptive reference files: their URLs no longer exist (HTTP 404), so the
#   tested image has an empty db/kaptive folder too (Kaptive is used only for
#   Klebsiella).
# config_file.txt is the tested image's file (it already points to the Kraken2
# database installed by the Dockerfile).
set -euo pipefail

ENV_DIR=/opt/conda/envs/tormes-1.3.0
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# CGE modules: repository, folder, commit
cge() {
    git clone -q "https://bitbucket.org/genomicepidemiology/$1" "$ENV_DIR/CGE-modules/$2"
    git -C "$ENV_DIR/CGE-modules/$2" checkout -q "$3"
    test "$(git -C "$ENV_DIR/CGE-modules/$2" rev-parse HEAD)" = "$3"
}
mkdir -p "$ENV_DIR/CGE-modules"
cge pointfinder.git        pointfinder                    b9d02190e33d63cbca04ad86599376fed57ba137
cge pointfinder_db.git     pointfinder/pointfinder_db     44ce624a806c6d2b70f7e39841a5f9cb4d9010aa
cge fimtyper.git           fimtyper                       311023be5c2caf8ad0673b15169ee02af5c46954
cge fimtyper_db            fimtyper/fimtyper_db           1dd6cdd40bdb967cdbcb049e6e835c1df840f4af
cge serotypefinder.git     serotypefinder                 d23ee1bcbb2dcf56d9ba3b00f6dbd23a583ad669
cge serotypefinder_db.git  serotypefinder/serotypefinder_db d9be114411a6561e8b5c43db292737c7275195f5

# roary2svg.pl from a fixed Roary commit
wget -q -O "$ENV_DIR/bin/roary2svg.pl" \
    https://raw.githubusercontent.com/sanger-pathogens/Roary/5cc13c8f26c5257204d7b88c7eb5e1db5031f9da/contrib/roary2svg/roary2svg.pl
echo "113e9cc542950c64e219f76e0ce93589  $ENV_DIR/bin/roary2svg.pl" | md5sum -c --quiet
chmod 755 "$ENV_DIR/bin/roary2svg.pl"

mkdir -p "$ENV_DIR/db/kaptive"

# R packages that tormes-setup adds with BiocManager (Bioconductor 3.10):
# ggtree 2.0.0 (conda) -> 2.0.4, and BiocVersion 3.10.1
mkdir -p /tmp/rpkgs
for p in ggtree_2.0.4 BiocVersion_3.10.1; do
    wget -q -P /tmp/rpkgs "https://bioconductor.org/packages/3.10/bioc/src/contrib/$p.tar.gz"
done
conda run -n tormes-1.3.0 R CMD INSTALL --library="$ENV_DIR/lib/R/library" \
    /tmp/rpkgs/BiocVersion_3.10.1.tar.gz /tmp/rpkgs/ggtree_2.0.4.tar.gz
rm -rf /tmp/rpkgs
conda run -n tormes-1.3.0 Rscript -e 'stopifnot(packageVersion("ggtree") == "2.0.4", packageVersion("BiocVersion") == "3.10.1")'

mkdir -p "$ENV_DIR/files"
cp "$HERE/config_file.txt" "$ENV_DIR/files/config_file.txt"

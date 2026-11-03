#!/usr/bin/env bash
# Download all FASTQ files (and run/sample metadata) for an SRA/ENA project.
# Uses the ENA mirror of SRA: files are already gzipped FASTQ, no sra-tools needed.
#
# Usage: ./download_PRJNA1198075.sh [PROJECT] [OUTDIR] [PARALLEL_JOBS]
set -euo pipefail

PROJECT="${1:-PRJNA1198075}"
OUTDIR="${2:-${PROJECT}_data}"
JOBS="${3:-4}"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

FIELDS="run_accession,sample_accession,experiment_accession,sample_alias,sample_title,library_layout,library_strategy,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
echo "Fetching metadata for $PROJECT ..."
curl -sSf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=${FIELDS}&format=tsv&download=true" \
    -o metadata.tsv

# Build a "url<TAB>md5" list, one line per FASTQ file
awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) col[$i]=i; next }
    { n=split($col["fastq_ftp"], u, ";"); split($col["fastq_md5"], m, ";");
      for (i=1;i<=n;i++) print "https://" u[i] "\t" m[i] }' metadata.tsv > files.tsv

# Sample groups: ENA leaves sample_alias/title empty, so fetch the BioSample
# records from NCBI. Sample names D_GP_xx = DOP (= PDO), ND_GP_xx = non-DOP (NONPDO).
echo "Fetching BioSample records and writing samplesheet.tsv ..."
SAMPLE_IDS="$(tail -n +2 metadata.tsv | cut -f2 | paste -sd,)"
curl -sSf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=${SAMPLE_IDS}&retmode=xml" \
    -o biosample.xml
# fastq_1/fastq_2 are relative to this directory, so the sheet still works when
# the directory is mounted somewhere else (e.g. inside a container).
python3 -I - biosample.xml metadata.tsv > samplesheet.tsv <<'EOF'
import sys, csv, xml.etree.ElementTree as ET
xml_file, meta_file = sys.argv[1:]
biosamples = {}
for b in ET.parse(xml_file).getroot().findall('BioSample'):
    name = next(i.text for i in b.findall('Ids/Id') if i.get('db_label') == 'Sample name')
    biosamples[b.get('accession')] = (name, b.findtext('Description/Title') or '')
rows = []
for r in csv.DictReader(open(meta_file), delimiter='\t'):
    name, title = biosamples[r['sample_accession']]
    group = 'NONPDO' if name.startswith('ND_') else 'PDO' if name.startswith('D_') else 'UNKNOWN'
    run = r['run_accession']
    rows.append((group != 'PDO', int(name.rsplit('_', 1)[1]),
                 [name, run, r['sample_accession'], title, group,
                  f'{run}_1.fastq.gz', f'{run}_2.fastq.gz', r['read_count']]))
print('\t'.join(['sample_name', 'run_accession', 'biosample', 'title', 'group',
                 'fastq_1', 'fastq_2', 'read_pairs']))
for *_, row in sorted(rows):
    print('\t'.join(row))
EOF
echo "Groups: $(tail -n +2 samplesheet.tsv | cut -f5 | sort | uniq -c | xargs)"

echo "$(tail -n +2 metadata.tsv | wc -l) runs, $(wc -l < files.tsv) files to download"

download_one() {
    url="$1"; md5="$2"; f="$(basename "$url")"
    if [[ -s "$f" ]] && echo "$md5  $f" | md5sum -c --status; then
        echo "OK (already present) $f"; return 0
    fi
    for attempt in 1 2 3 4 5 6; do
        # -C - resumes a partial download. A dropped connection keeps the
        # partial file so the next attempt resumes; a complete file with a
        # wrong checksum is corrupt and is downloaded again from scratch.
        if curl -sSfL --retry 10 -C - -o "$f" "$url"; then
            if echo "$md5  $f" | md5sum -c --status; then
                echo "OK $f"; return 0
            fi
            echo "Checksum mismatch for $f, downloading again" >&2; rm -f "$f"
        else
            echo "Connection dropped for $f (attempt $attempt), resuming" >&2; sleep 30
        fi
    done
    echo "FAILED $f" >&2; return 1
}
export -f download_one

# Download in parallel; xargs exits non-zero if any file failed
tr '\t' ' ' < files.tsv | xargs -P "$JOBS" -n 2 bash -c 'download_one "$0" "$1"'

echo "All files downloaded and verified in $(pwd)"

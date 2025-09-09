#!/bin/bash
set -Euo pipefail
shopt -s nullglob

# ===== Env =====
source /opt/conda/etc/profile.d/conda.sh
conda activate bioenv

# ===== Directories and parameters =====
input_dir="@00-01_dir1@"
output_dir00="@00-01_dir2@"   # non usato qui, ma lo mantengo per coerenza con la pipeline
output_dir01="@00-01_dir3@"
threads="@00-01_par1@"

mkdir -p "$output_dir00" "$output_dir01"

# ===== Log (nome corto) =====
ts="$(date +'%F')"   # YYYY-MM-DD
LOG="$output_dir01/00-01_adapterremoval_$ts.log"

log(){ echo "[$(date +'%F %T')] $*" | tee -a "$LOG"; }
fail_msg(){ echo "[$(date +'%F %T')] ERROR: $*" | tee -a "$LOG" >&2; }

log "== Starting step 00-01 (AdapterRemoval) =="
log "Input: $input_dir"
log "Output cleaned: $output_dir01"
log "Threads: $threads"
echo >> "$LOG"

# ===== Results collections =====
ok_samples=()
failed_missing_pair=()
failed_integrity=()
failed_adapterremoval=()

log "Starting AdapterRemoval..."

# Loop: gestisce sia *_1.fq.gz che *_L1.fq.gz ed evita *_cleaned_*
for file_L1 in "$input_dir"/*_1.fq.gz "$input_dir"/*_L1.fq.gz; do
  [[ -e "$file_L1" ]] || continue

  base="$(basename "$file_L1")"
  # skip già puliti
  if [[ "$base" == *_cleaned_* ]]; then
    log "Skip already-cleaned file: $base"
    continue
  fi

  # ricava il prefisso sample e normalizza la ricerca del mate
  if [[ "$base" == *_1.fq.gz ]]; then
    file_base="${file_L1%_1.fq.gz}"
    file_L2="${file_base}_2.fq.gz"
  else
    file_base="${file_L1%_L1.fq.gz}"
    file_L2="${file_base}_L2.fq.gz"
  fi

  sample_name="$(basename "$file_base")"

  # 1) Check pair presente
  if [[ ! -f "$file_L2" ]]; then
    # prova l'altra convenzione (se era _1, cerca _L2; se era _L1, cerca _2)
    alt_L2="${file_base}_2.fq.gz"
    alt2_L2="${file_base}_L2.fq.gz"
    if [[ -f "$alt_L2" ]]; then
      file_L2="$alt_L2"
    elif [[ -f "$alt2_L2" ]]; then
      file_L2="$alt2_L2"
    else
      fail_msg "Missing paired file for $base (tried: $file_L2, $alt_L2, $alt2_L2)"
      failed_missing_pair+=("$sample_name")
      continue
    fi
  fi

  # 2) Integrità gzip
  if ! gzip -t "$file_L1" 2>>"$LOG"; then
    fail_msg "Integrity check failed (gzip -t): $(basename "$file_L1")"
    failed_integrity+=("$sample_name")
    continue
  fi
  if ! gzip -t "$file_L2" 2>>"$LOG"; then
    fail_msg "Integrity check failed (gzip -t): $(basename "$file_L2")"
    failed_integrity+=("$sample_name")
    continue
  fi

  # 3) AdapterRemoval (output nella cartella di output)
  out1="$output_dir01/${sample_name}_cleaned_L1.fq.gz"
  out2="$output_dir01/${sample_name}_cleaned_L2.fq.gz"
  ar_base="$output_dir01/${sample_name}"

  log "--> Sample $sample_name: running AdapterRemoval"
  if AdapterRemoval \
      --file1 "$file_L1" --file2 "$file_L2" \
      --output1 "$out1" \
      --output2 "$out2" \
      --basename "$ar_base" \
      --discarded "$output_dir01/${sample_name}_discarded.fq.gz" \
      --settings "$output_dir01/${sample_name}.settings" \
      --threads "$threads" \
      --gzip \
      --minlength @00-01_par2@ \
      --trimqualities \
      --minquality @00-01_par3@ \
      --trimns \
      --maxns @00-01_par4@ \
      --trim5p @00-01_par5@ \
      --trim3p @00-01_par6@ >>"$LOG" 2>&1; then
    log "OK AdapterRemoval: $sample_name"
    ok_samples+=("$sample_name")
  else
    fail_msg "AdapterRemoval failed for $sample_name (see log)"
    failed_adapterremoval+=("$sample_name")
    # continua coi successivi
  fi
done

# ===== Summary =====
total=$(( ${#ok_samples[@]} + ${#failed_missing_pair[@]} + ${#failed_integrity[@]} + ${#failed_adapterremoval[@]} ))
log ""
log "== Summary =="
log "Total samples processed: $total"
log "   OK: ${#ok_samples[@]}"
log "   Missing pair: ${#failed_missing_pair[@]}"
log "   Gzip integrity failed: ${#failed_integrity[@]}"
log "   AdapterRemoval failed: ${#failed_adapterremoval[@]}"

print_list(){
  local title="$1"; shift
  local -a arr=("$@")
  if (( ${#arr[@]} > 0 )); then
    log "$title (${#arr[@]}):"
    for s in "${arr[@]}"; do log "  - $s"; done
  else
    log "$title: none"
  fi
}

print_list "OK" "${ok_samples[@]}"
print_list "MISSING PAIR" "${failed_missing_pair[@]}"
print_list "FAILED INTEGRITY" "${failed_integrity[@]}"
print_list "AR FAILED" "${failed_adapterremoval[@]}"

echo
echo "=== RESULT: ${#ok_samples[@]} out of $total OK ==="
echo "Detailed log: $LOG"

# Exit code: 0 se tutti OK, 1 se c'è almeno un fallimento
if (( ${#failed_missing_pair[@]} + ${#failed_integrity[@]} + ${#failed_adapterremoval[@]} > 0 )); then
  exit 1
fi
exit 0

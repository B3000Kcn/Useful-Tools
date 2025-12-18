#!/usr/bin/env bash
# autocycler.sh (stable, parameterized)
# Usage:
#   ./autocycler.sh reads.fq -t 64 -j 1 -i assemblies_dir -o out_dir --max_contigs 30
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: autocycler.sh <reads.fq> [OPTIONS]

Options:
  -t, --threads <int>      Threads used ONLY for autocycler compress (default: 8)
  -j, --jobs <int>         Parallel cluster processing jobs for trim+resolve (default: 1)
  -i, --input_dir <path>   Directory containing input assemblies (*.fa/*.fasta). If set, skip generation.
  -o, --out_dir <path>     Output directory (default: autocycler_out)
  --max_contigs <int>      Max contigs threshold (default: 25). Also used for pre-filter.
  -h, --help               Show help
EOF
}

if [[ $# -lt 1 ]]; then usage; exit 1; fi
reads=$(realpath "$1")
shift

threads=8
jobs=1
out_dir="autocycler_out"
max_contigs=25
input_assemblies=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--threads) threads="$2"; shift 2 ;;
    -j|--jobs) jobs="$2"; shift 2 ;;
    -i|--input_dir) input_assemblies=$(realpath "$2"); shift 2 ;;
    -o|--out_dir) out_dir=$(realpath "$2"); shift 2 ;;
    --max_contigs) max_contigs="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

mkdir -p "$out_dir"
log_file="$out_dir/autocycler.stderr"

echo "=========================================="
echo " Autocycler Wrapper (Stable Param Version)"
echo "=========================================="
echo "Reads:       $reads"
echo "Out:         $out_dir"
echo "Threads:     $threads  (only for compress)"
echo "Jobs:        $jobs     (parallel clusters trim+resolve)"
echo "Max contigs: $max_contigs"
echo "=========================================="

# ---------- Mode: Consumer (external assemblies) ----------
if [[ -n "$input_assemblies" ]]; then
  assemblies_dir="$out_dir/assemblies_work"
  mkdir -p "$assemblies_dir"
  echo "[Mode] Consumer. Copy assemblies from: $input_assemblies"

  shopt -s nullglob
  cp "$input_assemblies"/*.fasta "$assemblies_dir"/ 2>/dev/null || true
  cp "$input_assemblies"/*.fa    "$assemblies_dir"/ 2>/dev/null || true
  shopt -u nullglob

  # Add weight to all copied assemblies (safe: working dir)
  echo "-> Applying consensus weight=2 to assemblies..."
  shopt -s nullglob
  for f in "$assemblies_dir"/*.{fasta,fa}; do
    [[ -e "$f" ]] || continue
    sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
  done
  shopt -u nullglob
else
  echo "❌ This script version expects -i/--input_dir (consumer mode)."
  echo "   If you want producer mode (subsample+assemble), tell me and I’ll merge it back."
  exit 1
fi

# ---------- Pre-filter (The Sieve) ----------
echo -e "\n-> [Pre-check] Filtering assemblies with contigs > $max_contigs ..."
mkdir -p "$out_dir/discarded_fragmented"
removed=0

shopt -s nullglob
for f in "$assemblies_dir"/*.{fasta,fa}; do
  [[ -s "$f" ]] || continue
  n=$(grep -c "^>" "$f" || true)
  if (( n > max_contigs )); then
    echo "   [Discard] $(basename "$f") has $n contigs"
    mv "$f" "$out_dir/discarded_fragmented/"
    ((++removed))
  fi
done
shopt -u nullglob

echo "-> Discarded: $removed"

if [[ "$(ls -A "$assemblies_dir" 2>/dev/null | wc -l)" -eq 0 ]]; then
  echo "❌ Error: all assemblies discarded. Increase --max_contigs or fix assemblies."
  exit 1
fi

# ---------- Step 3: compress (threads supported) ----------
echo -e "\n-> Step 3: autocycler compress (threads=$threads)"
autocycler compress \
  --assemblies_dir "$assemblies_dir" \
  --autocycler_dir "$out_dir" \
  --max_contigs "$max_contigs" \
  --threads "$threads" \
  2>>"$log_file"

# ---------- Step 4: cluster (NO --threads in your version) ----------
echo -e "\n-> Step 4: autocycler cluster (max_contigs=$max_contigs)"
autocycler cluster \
  --autocycler_dir "$out_dir" \
  --max_contigs "$max_contigs" \
  2>>"$log_file"

# ---------- Step 5+6: trim + resolve (NO --threads for resolve in your version) ----------
echo -e "\n-> Step 5+6: trim + resolve (jobs=$jobs)"
shopt -s nullglob
clusters=( "$out_dir"/clustering/qc_pass/cluster_* )
shopt -u nullglob

if [[ ${#clusters[@]} -eq 0 ]]; then
  echo "❌ Error: No qc_pass clusters. Check $log_file"
  exit 1
fi

process_one_cluster() {
  local c="$1"
  echo "   Processing $(basename "$c")"
  # Use --cluster_dir to match your resolve usage output
  autocycler trim   --cluster_dir "$c" 2>>"$log_file"
  autocycler resolve --cluster_dir "$c" 2>>"$log_file"
}

# parallelize across clusters if jobs > 1
if (( jobs > 1 )) && (( ${#clusters[@]} > 1 )); then
  # simple bash job control, no GNU parallel dependency
  running=0
  for c in "${clusters[@]}"; do
    process_one_cluster "$c" &
    ((running++))
    if (( running >= jobs )); then
      wait -n
      ((running--))
    fi
  done
  wait
else
  for c in "${clusters[@]}"; do
    process_one_cluster "$c"
  done
fi

# ---------- Step 7: combine ----------
echo -e "\n-> Step 7: autocycler combine"
autocycler combine \
  -a "$out_dir" \
  -i "$out_dir"/clustering/qc_pass/cluster_*/5_final.gfa \
  2>>"$log_file"

echo -e "\n✅ Done. Final: $out_dir/assembly.fasta"

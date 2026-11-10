#!/usr/bin/env bash
# Plastizyme screen of MAGs (Sci Data revision: YJ vs JJ field soils).
# Mirrors modules/prodigal.nf, plastizyme_qc.nf and diamond.nf of this repository.
# Usage: run_plastizyme_screen.sh <bins_root> <refs_dir> <outdir> [threads]
#   bins_root/<SAMPLE>/*.fasta  (one sub-folder per sample, e.g. YJ, JJ)
#   refs_dir: PLA_pazy_proteins.fasta PLA_PlasticDB.fasta PBAT_pazy_proteins.fasta PBAT_PlasticDB.fasta
set -euo pipefail
BINS=$1; REFS=$2; OUT=$3; T=${4:-4}
PRODIGAL=${PRODIGAL:-prodigal}; DIAMOND=${DIAMOND:-diamond}; CDHIT=${CDHIT:-cd-hit}
mkdir -p "$OUT"/{prodigal,qc,diamond}

# 1. Gene calling (single mode, meta fallback), headers prefixed SAMPLE__BIN__
: > "$OUT/prodigal/all_bins.faa"
for sdir in "$BINS"/*/; do
  s=$(basename "$sdir")
  for f in "$sdir"*.fa "$sdir"*.fasta; do
    [ -e "$f" ] || continue
    b=$(basename "$f"); b=${b%.*}; id="${s}__${b}"
    if ! $PRODIGAL -i "$f" -a "$OUT/prodigal/$id.faa" -d "$OUT/prodigal/$id.ffn" \
         -o "$OUT/prodigal/$id.gff" -f gff -p single -q 2> "$OUT/prodigal/$id.log"; then
      $PRODIGAL -i "$f" -a "$OUT/prodigal/$id.faa" -d "$OUT/prodigal/$id.ffn" \
         -o "$OUT/prodigal/$id.gff" -f gff -p meta -q
    fi
    sed -i "s/^>/>${id}__/" "$OUT/prodigal/$id.faa" "$OUT/prodigal/$id.ffn"
    cat "$OUT/prodigal/$id.faa" >> "$OUT/prodigal/all_bins.faa"
  done
done
echo "Proteins: $(grep -c '^>' "$OUT/prodigal/all_bins.faa")"

# 2. Reference QC: merge PAZy + PlasticDB, keep >=100 aa, CD-HIT 95 %, DIAMOND db
for db in PLA PBAT; do
  q="$OUT/qc/$db"
  cat "$REFS/${db}_pazy_proteins.fasta" "$REFS/${db}_PlasticDB.fasta" > "$q.merged.faa"
  awk 'BEGIN{seq="";h=""} /^>/{if(h!="" && length(seq)>=100) print h"\n"seq; h=$0; seq=""; next} {seq=seq$0}
       END{if(h!="" && length(seq)>=100) print h"\n"seq}' "$q.merged.faa" > "$q.filtered.faa"
  $CDHIT -i "$q.filtered.faa" -o "$q.clean.faa" -c 0.95 -n 5 -T "$T" -M 2000 > "$q.cdhit.log"
  {
    echo "$db merged: $(grep -c '^>' "$q.merged.faa")"
    echo "$db >=100 aa: $(grep -c '^>' "$q.filtered.faa")"
    echo "$db after CD-HIT 95%: $(grep -c '^>' "$q.clean.faa")"
  } | tee "$q.qc_report.txt"
  $DIAMOND makedb --in "$q.clean.faa" --db "$q" --threads "$T" --quiet

# 3. DIAMOND blastp (same thresholds as modules/diamond.nf)
  $DIAMOND blastp --query "$OUT/prodigal/all_bins.faa" --db "$q.dmnd" \
    --out "$OUT/diamond/$db.diamond.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
    --evalue 1e-5 --id 30 --query-cover 70 --threads "$T" --more-sensitive --quiet
  echo "$db hits: $(wc -l < "$OUT/diamond/$db.diamond.tsv")"
done

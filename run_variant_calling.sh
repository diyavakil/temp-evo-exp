#!/usr/bin/env bash
# run_variant_calling.sh
#
# Local re-creation of the small-variant (and optional structural-variant) pipeline that Plasmidsaurus
# describes for its 'Call Small Variants' step:
#   1. Filtlong      reduce raw ONT reads to at most ~100x coverage
#   2. minimap2      align the filtered reads to the reference
#   3. Clair3        call SNPs and small indels (bacteria-tuned model)
#   4. bcftools      filter low-confidence and low-depth calls
#   5. SnpEff        annotate and predict effects (needs the reference as GenBank/GBFF)
#   6. (optional) Sniffles2 or assembly-to-assembly alignment for larger structural variants
#   7. vcf_to_table.py   write per-replicon CSV tables that subtract_variants.py can read
#
# *** STATUS: UNTESTED RECONSTRUCTION. ***
# Written from the published description of the Plasmidsaurus pipeline. The exact Filtlong options, Clair3
# model, Clair3 options and filter thresholds they use are not stated, so the defaults below are placeholders.
# Before relying on it, run it on a sample you already have Plasmidsaurus calls for (for example 18C ESL1) and
# compare the two call sets (see docs/VARIANT_CALLING_LOCAL.md).
#
# USAGE
#   scripts/run_variant_calling.sh -s 18C_ESL1 -r ref.fasta -g ref.gbff -q reads.fastq.gz \
#       -m /path/to/clair3_bacteria_model -o results/18C_ESL1 [-t 8] [-c 100] [-Q 10] [-D 20] \
#       [-M NC_006840.2=chr1,NC_006841.2=chr2,NC_006842.1=pES100] [-V sniffles|asm|none] [-a assembly.fasta]
#
#   -s  sample name (used in file names)                      -r  reference FASTA
#   -g  reference GenBank/GBFF with annotation (for SnpEff)   -q  raw reads (FASTQ or FASTQ.gz)
#   -m  Clair3 model directory (bacteria-tuned, matching your basecalling model)
#   -o  output directory   -t threads (8)   -c max coverage (100)   -Q min QUAL (10)   -D min depth (20)
#   -M  contig=replicon map for the CSV tables (chr1/chr2/pES100 names)
#   -V  structural variants: sniffles (default), asm (assembly vs reference), none      -a assembly FASTA (asm mode)
#
# Optional environment variables: CLAIR3_EXTRA (extra Clair3 options, e.g. "--haploid_precise"),
# SNPEFF (path to the snpEff command, default 'snpEff').

set -euo pipefail

THREADS=8; MAXCOV=100; MIN_QUAL=10; MIN_DP=20; SVMODE=sniffles
REPMAP="NC_006840.2=chr1,NC_006841.2=chr2,NC_006842.1=pES100"
SNPEFF=${SNPEFF:-snpEff}; CLAIR3_EXTRA=${CLAIR3_EXTRA:-}
SAMPLE=""; REF=""; GBFF=""; READS=""; MODEL=""; OUT=""; ASM=""

usage() { sed -n '2,34p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts "s:r:g:q:m:o:t:c:Q:D:M:V:a:h" opt; do
  case $opt in
    s) SAMPLE=$OPTARG;; r) REF=$OPTARG;; g) GBFF=$OPTARG;; q) READS=$OPTARG;; m) MODEL=$OPTARG;;
    o) OUT=$OPTARG;; t) THREADS=$OPTARG;; c) MAXCOV=$OPTARG;; Q) MIN_QUAL=$OPTARG;; D) MIN_DP=$OPTARG;;
    M) REPMAP=$OPTARG;; V) SVMODE=$OPTARG;; a) ASM=$OPTARG;; h|*) usage;;
  esac
done
for v in SAMPLE REF READS MODEL OUT; do [[ -n "${!v}" ]] || { echo "missing required option for $v" >&2; usage; }; done
[[ -f "$REF" && -f "$READS" ]] || { echo "reference or reads file not found" >&2; exit 1; }
for tool in filtlong minimap2 samtools bcftools run_clair3.sh; do
  command -v "$tool" >/dev/null || { echo "required tool not found on PATH: $tool" >&2; exit 1; }
done

mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
LOG="$OUT/${SAMPLE}.versions.txt"
{ echo "# tool versions ($(date -u +%Y-%m-%dT%H:%MZ))"
  filtlong --version 2>&1 | head -1; minimap2 --version; samtools --version | head -1; bcftools --version | head -1
  command -v sniffles >/dev/null && sniffles --version 2>&1 | head -1 || true
  echo "MAXCOV=$MAXCOV MIN_QUAL=$MIN_QUAL MIN_DP=$MIN_DP SVMODE=$SVMODE CLAIR3_EXTRA='$CLAIR3_EXTRA'"; } > "$LOG"

# ---- 1. quality filtering to at most MAXCOV x coverage -------------------------------------------------
GENOME=$(grep -v '>' "$REF" | tr -d '\n\r ' | wc -c)
TARGET=$(( GENOME * MAXCOV ))
echo "[1/7] Filtlong: genome ${GENOME} bp, keeping best reads up to ${TARGET} bases"
filtlong --target_bases "$TARGET" "$READS" | gzip > "$OUT/${SAMPLE}.filtered.fastq.gz"

# ---- 2. alignment -----------------------------------------------------------------------------------------
echo "[2/7] minimap2 alignment"
samtools faidx "$REF"
minimap2 -ax map-ont -t "$THREADS" "$REF" "$OUT/${SAMPLE}.filtered.fastq.gz" \
  | samtools sort -@ "$THREADS" -o "$OUT/${SAMPLE}.bam" -
samtools index "$OUT/${SAMPLE}.bam"

# ---- 3. variant calling -----------------------------------------------------------------------------------
echo "[3/7] Clair3"
run_clair3.sh --bam_fn="$OUT/${SAMPLE}.bam" --ref_fn="$REF" --threads="$THREADS" --platform=ont \
  --model_path="$MODEL" --output="$OUT/clair3" --include_all_ctgs $CLAIR3_EXTRA

# ---- 4. filtering ---------------------------------------------------------------------------------------
echo "[4/7] bcftools filter (QUAL >= ${MIN_QUAL}, depth >= ${MIN_DP}, PASS)"
bcftools view -f PASS -i "QUAL>=${MIN_QUAL} && FMT/DP>=${MIN_DP}" "$OUT/clair3/merge_output.vcf.gz" \
  -Oz -o "$OUT/${SAMPLE}.filtered.vcf.gz"
bcftools index -t "$OUT/${SAMPLE}.filtered.vcf.gz"

# ---- 5. annotation (SnpEff, database built from the reference GenBank file) --------------------------------
ANNOTATED="$OUT/${SAMPLE}.filtered.vcf.gz"
if [[ -n "$GBFF" ]]; then
  echo "[5/7] SnpEff"
  command -v "$SNPEFF" >/dev/null || { echo "snpEff not found (set SNPEFF=...)" >&2; exit 1; }
  GID="${SAMPLE}_ref"; SDIR="$OUT/snpeff"; mkdir -p "$SDIR/data/$GID"
  cp "$GBFF" "$SDIR/data/$GID/genes.gbk"
  printf 'data.dir = %s/\n%s.genome : %s\n' "$SDIR/data" "$GID" "$GID" > "$SDIR/snpEff.config"
  "$SNPEFF" build -genbank -v -c "$SDIR/snpEff.config" "$GID"
  "$SNPEFF" -c "$SDIR/snpEff.config" "$GID" "$OUT/${SAMPLE}.filtered.vcf.gz" > "$OUT/${SAMPLE}.annotated.vcf"
  ANNOTATED="$OUT/${SAMPLE}.annotated.vcf"
else
  echo "[5/7] no GenBank reference given (-g): skipping SnpEff; tables will have no gene annotation"
fi

# ---- 6. structural variants (optional) -----------------------------------------------------------------------
case "$SVMODE" in
  sniffles)
    echo "[6/7] Sniffles2"
    command -v sniffles >/dev/null || { echo "sniffles not found (use -V none to skip)" >&2; exit 1; }
    sniffles --input "$OUT/${SAMPLE}.bam" --reference "$REF" --threads "$THREADS" --vcf "$OUT/${SAMPLE}.sv.vcf.gz" ;;
  asm)
    echo "[6/7] assembly-to-assembly alignment (minimap2 + paftools.js)"
    [[ -f "$ASM" ]] || { echo "-a assembly FASTA required for -V asm" >&2; exit 1; }
    minimap2 -cx asm5 -t "$THREADS" --cs "$REF" "$ASM" | sort -k6,6 -k8,8n \
      | paftools.js call -f "$REF" -l 50 -L 50 - > "$OUT/${SAMPLE}.asm_calls.txt"
    echo "    raw paftools output written; keep calls >= 50 bp when you read it" ;;
  none) echo "[6/7] structural variants skipped" ;;
  *) echo "unknown -V mode: $SVMODE" >&2; exit 1 ;;
esac

# ---- 7. tables ---------------------------------------------------------------------------------------------
echo "[7/7] writing CSV tables"
python3 "$HERE/vcf_to_table.py" "$ANNOTATED" --sample "$SAMPLE" --out-dir "$OUT/tables" --rep-map "$REPMAP" --pass-only
echo "done. Small-variant tables: $OUT/tables ; tool versions: $LOG"

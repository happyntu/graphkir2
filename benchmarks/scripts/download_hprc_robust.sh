#!/usr/bin/env bash
# Robust HPRC FASTQ downloader with per-sample error isolation,
# pigz compression, skip-if-exists, and automatic MD5 recording.

set -uo pipefail

FASTQ_ROOT=/mypool/KIR_graph/data/hprc_fastq
SRA_CACHE=/mypool/KIR_graph/tmp/sra
TMP_DIR=/mypool/KIR_graph/tmp/fasterq
LOG_DIR=/mypool/KIR_graph/logs
MD5_MANIFEST=/mypool/KIR_graph/data/hprc_fastq/md5_manifest.tsv
THREADS=8
PIGZ_THREADS=4
DRY_RUN=false
SINGLE_SAMPLE=""

export PATH=/mypool/KIR_graph/envs/graphkir_env/bin:$PATH

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=true; shift ;;
        --sample) SINGLE_SAMPLE="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

mkdir -p "$FASTQ_ROOT" "$SRA_CACHE" "$TMP_DIR" "$LOG_DIR"

# Initialize MD5 manifest if it does not exist
if [[ ! -f "$MD5_MANIFEST" ]]; then
    printf 'sample_id\taccession\tread\tfile\tmd5\tbytes\ttimestamp\n' > "$MD5_MANIFEST"
fi

PLAN_TSV=/mypool/KIR_graph/tmp/download-plan/hprc_fastq_download_plan.tsv

record_md5() {
    local sample_id="$1"
    local accession="$2"
    local read_label="$3"
    local filepath="$4"
    local log="$5"

    if [[ ! -f "$filepath" ]]; then
        echo "  WARN: $filepath not found, skipping MD5" | tee -a "$log"
        return
    fi

    local md5_line
    md5_line=$(md5sum "$filepath")
    local md5_val="${md5_line%% *}"
    local bytes
    bytes=$(stat -c%s "$filepath")
    local ts
    ts=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_id" "$accession" "$read_label" \
        "$(basename "$filepath")" "$md5_val" "$bytes" "$ts" \
        >> "$MD5_MANIFEST"

    echo "  MD5 ${read_label}: ${md5_val} [${bytes} bytes]" | tee -a "$log"
}

download_sample() {
    local sample_id="$1"
    local accession="$2"
    local r1="$FASTQ_ROOT/${accession}_1.fastq.gz"
    local r2="$FASTQ_ROOT/${accession}_2.fastq.gz"
    local sra_path="$SRA_CACHE/$accession/$accession.sra"
    local sample_tmp="$TMP_DIR/$accession"
    local log="$LOG_DIR/download_${accession}.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $sample_id / $accession" | tee -a "$log"

    # Skip if both final FASTQs exist and are non-empty
    if [[ -s "$r1" && -s "$r2" ]]; then
        echo "  SKIP: both reads exist" | tee -a "$log"
        return 0
    fi

    if $DRY_RUN; then
        echo "  DRY-RUN: would download $accession" | tee -a "$log"
        return 0
    fi

    mkdir -p "$sample_tmp"

    # Prefetch if SRA not cached
    if [[ ! -f "$sra_path" ]]; then
        echo "  prefetch $accession ..." | tee -a "$log"
        prefetch "$accession" --check-rs no --max-size 100G \
            --output-directory "$SRA_CACHE" 2>&1 | tee -a "$log"
    else
        echo "  SRA cached: $sra_path" | tee -a "$log"
    fi

    # Verify SRA file exists after prefetch
    if [[ ! -f "$sra_path" ]]; then
        echo "  ERROR: SRA file not found after prefetch: $sra_path" | tee -a "$log"
        return 1
    fi

    # fasterq-dump
    echo "  fasterq-dump $accession ..." | tee -a "$log"
    fasterq-dump "$sra_path" --split-files --threads "$THREADS" \
        --outdir "$sample_tmp" 2>&1 | tee -a "$log"

    # Compress with pigz and record MD5
    local raw1="$sample_tmp/${accession}_1.fastq"
    local raw2="$sample_tmp/${accession}_2.fastq"

    if [[ -f "$raw1" ]]; then
        echo "  pigz R1 ..." | tee -a "$log"
        pigz -p "$PIGZ_THREADS" -c "$raw1" > "$r1"
        rm -f "$raw1"
    fi

    if [[ -f "$raw2" ]]; then
        echo "  pigz R2 ..." | tee -a "$log"
        pigz -p "$PIGZ_THREADS" -c "$raw2" > "$r2"
        rm -f "$raw2"
    fi

    # Record MD5 for both reads
    echo "  Computing MD5 checksums ..." | tee -a "$log"
    record_md5 "$sample_id" "$accession" "R1" "$r1" "$log"
    record_md5 "$sample_id" "$accession" "R2" "$r2" "$log"

    # Clean up SRA cache for this sample to save disk
    echo "  Cleaning SRA cache for $accession ..." | tee -a "$log"
    rm -rf "$SRA_CACHE/$accession"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE $sample_id / $accession" | tee -a "$log"
}

# Parse plan TSV (skip header)
SUCCESS=0
FAIL=0

tail -n+2 "$PLAN_TSV" | while IFS=$'\t' read -r sid accession rest; do
    if [[ -n "$SINGLE_SAMPLE" && "$accession" != "$SINGLE_SAMPLE" ]]; then
        continue
    fi

    if download_sample "$sid" "$accession"; then
        ((SUCCESS++)) || true
    else
        echo "  FAILED: $sid / $accession" | tee -a "$LOG_DIR/download_${accession}.log"
        ((FAIL++)) || true
    fi
done

echo ""
echo "Download complete. Success: $SUCCESS, Failed: $FAIL"
echo "FASTQ root: $FASTQ_ROOT"
echo "MD5 manifest: $MD5_MANIFEST"
echo "Logs: $LOG_DIR/download_*.log"

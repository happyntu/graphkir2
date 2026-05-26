#!/usr/bin/env bash
# Parallel HPRC FASTQ downloader.
# Runs up to PARALLEL_JOBS concurrent sample pipelines (prefetch + fasterq-dump + pigz).
# Each sample is isolated: per-sample log, per-sample MD5 file, independent error handling.
# Safe to restart: skips samples where both final FASTQs already exist;
# prefetch resumes from .tmp files automatically.

set -uo pipefail

FASTQ_ROOT=/mypool/KIR_graph/data/hprc_fastq
SRA_CACHE=/mypool/KIR_graph/tmp/sra
TMP_DIR=/mypool/KIR_graph/tmp/fasterq
LOG_DIR=/mypool/KIR_graph/logs
MD5_DIR=/mypool/KIR_graph/logs/md5
MD5_MANIFEST=/mypool/KIR_graph/data/hprc_fastq/md5_manifest.tsv
THREADS=8
PIGZ_THREADS=4
PARALLEL_JOBS=3
DRY_RUN=false
SINGLE_SAMPLE=""

export PATH=/mypool/KIR_graph/envs/graphkir_env/bin:$PATH

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=true; shift ;;
        --sample) SINGLE_SAMPLE="$2"; shift 2 ;;
        --jobs|-j) PARALLEL_JOBS="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

mkdir -p "$FASTQ_ROOT" "$SRA_CACHE" "$TMP_DIR" "$LOG_DIR" "$MD5_DIR"

# Initialize MD5 manifest if it does not exist
if [[ ! -f "$MD5_MANIFEST" ]]; then
    printf 'sample_id\taccession\tread\tfile\tmd5\tbytes\ttimestamp\n' > "$MD5_MANIFEST"
fi

PLAN_TSV=/mypool/KIR_graph/tmp/download-plan/hprc_fastq_download_plan.tsv

# Per-thread count for fasterq-dump to avoid oversubscription
FQ_THREADS=$(( THREADS / PARALLEL_JOBS ))
[[ $FQ_THREADS -lt 2 ]] && FQ_THREADS=2

download_sample() {
    local sample_id="$1"
    local accession="$2"
    local r1="$FASTQ_ROOT/${accession}_1.fastq.gz"
    local r2="$FASTQ_ROOT/${accession}_2.fastq.gz"
    local sra_path="$SRA_CACHE/$accession/$accession.sra"
    local sample_tmp="$TMP_DIR/$accession"
    local log="$LOG_DIR/download_${accession}.log"
    local md5_file="$MD5_DIR/${accession}.md5.tsv"

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

    # Remove stale lock files left by interrupted prefetch
    local lock_file="$SRA_CACHE/$accession/$accession.sra.lock"
    if [[ -f "$lock_file" ]]; then
        echo "  Removing stale lock: $lock_file" | tee -a "$log"
        rm -f "$lock_file"
    fi

    # Prefetch if SRA not cached (resumes from .tmp automatically)
    if [[ ! -f "$sra_path" ]]; then
        echo "  prefetch $accession ..." | tee -a "$log"
        prefetch "$accession" --max-size 100G \
            --output-directory "$SRA_CACHE" 2>&1 | tee -a "$log"
    else
        echo "  SRA cached: $sra_path" | tee -a "$log"
    fi

    # Verify SRA file exists after prefetch
    if [[ ! -f "$sra_path" ]]; then
        echo "  ERROR: SRA file not found after prefetch: $sra_path" | tee -a "$log"
        return 1
    fi

    # fasterq-dump with per-job thread allocation
    echo "  fasterq-dump $accession (threads=$FQ_THREADS) ..." | tee -a "$log"
    fasterq-dump "$sra_path" --split-files --threads "$FQ_THREADS" \
        --outdir "$sample_tmp" 2>&1 | tee -a "$log"

    # Compress with pigz
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

    # Record MD5 to per-sample file (no locking needed)
    echo "  Computing MD5 checksums ..." | tee -a "$log"
    for read_label in R1 R2; do
        local filepath
        [[ "$read_label" == "R1" ]] && filepath="$r1" || filepath="$r2"

        if [[ ! -f "$filepath" ]]; then
            echo "  WARN: $filepath not found, skipping MD5" | tee -a "$log"
            continue
        fi

        local md5_val bytes ts
        md5_val=$(md5sum "$filepath" | cut -d' ' -f1)
        bytes=$(stat -c%s "$filepath")
        ts=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$sample_id" "$accession" "$read_label" \
            "$(basename "$filepath")" "$md5_val" "$bytes" "$ts" \
            >> "$md5_file"

        echo "  MD5 ${read_label}: ${md5_val} [${bytes} bytes]" | tee -a "$log"
    done

    # Clean up SRA cache for this sample to save disk
    echo "  Cleaning SRA cache for $accession ..." | tee -a "$log"
    rm -rf "$SRA_CACHE/$accession"

    # Clean up sample tmp dir
    rmdir "$sample_tmp" 2>/dev/null || true

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE $sample_id / $accession" | tee -a "$log"
}

export -f download_sample
export FASTQ_ROOT SRA_CACHE TMP_DIR LOG_DIR MD5_DIR THREADS PIGZ_THREADS FQ_THREADS DRY_RUN

echo "=== HPRC Parallel Download ==="
echo "  Jobs: $PARALLEL_JOBS"
echo "  Threads/job (fasterq-dump): $FQ_THREADS"
echo "  Plan: $PLAN_TSV"
echo "  FASTQ root: $FASTQ_ROOT"
echo ""

# Build job list from plan TSV (skip header)
JOB_LIST=$(mktemp)
tail -n+2 "$PLAN_TSV" | while IFS=$'\t' read -r sid accession rest; do
    if [[ -n "$SINGLE_SAMPLE" && "$accession" != "$SINGLE_SAMPLE" ]]; then
        continue
    fi
    echo "$sid $accession"
done > "$JOB_LIST"

TOTAL=$(wc -l < "$JOB_LIST")
echo "Total samples to process: $TOTAL"
echo ""

# Run jobs with controlled parallelism
RUNNING=0
SUCCESS=0
FAIL=0
declare -A PIDS

while IFS=' ' read -r sid accession; do
    # Wait if we've hit the parallel limit
    while [[ $RUNNING -ge $PARALLEL_JOBS ]]; do
        # Wait for any child to finish
        wait -n 2>/dev/null || true
        # Reap finished jobs
        for pid in "${!PIDS[@]}"; do
            if ! kill -0 "$pid" 2>/dev/null; then
                wait "$pid" 2>/dev/null
                local_exit=$?
                if [[ $local_exit -eq 0 ]]; then
                    ((SUCCESS++)) || true
                else
                    ((FAIL++)) || true
                    echo "  FAILED: ${PIDS[$pid]}" | tee -a "$LOG_DIR/parallel_summary.log"
                fi
                unset "PIDS[$pid]"
                ((RUNNING--)) || true
            fi
        done
    done

    # Launch background job
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] LAUNCH $sid / $accession (slot $((RUNNING+1))/$PARALLEL_JOBS)"
    download_sample "$sid" "$accession" &
    PIDS[$!]="$sid/$accession"
    ((RUNNING++)) || true

done < "$JOB_LIST"

# Wait for remaining jobs
for pid in "${!PIDS[@]}"; do
    wait "$pid" 2>/dev/null
    local_exit=$?
    if [[ $local_exit -eq 0 ]]; then
        ((SUCCESS++)) || true
    else
        ((FAIL++)) || true
        echo "  FAILED: ${PIDS[$pid]}" | tee -a "$LOG_DIR/parallel_summary.log"
    fi
done

rm -f "$JOB_LIST"

# Merge per-sample MD5 files into manifest
echo ""
echo "Merging MD5 records ..."
for f in "$MD5_DIR"/*.md5.tsv; do
    [[ -f "$f" ]] && cat "$f" >> "$MD5_MANIFEST"
done

echo ""
echo "=== Download Complete ==="
echo "  Success: $SUCCESS"
echo "  Failed:  $FAIL"
echo "  FASTQ root: $FASTQ_ROOT"
echo "  MD5 manifest: $MD5_MANIFEST"
echo "  Logs: $LOG_DIR/download_*.log"

#!/usr/bin/env bash
# Batch Geny runner for 44 HPRC samples.
# Runs PARALLEL_JOBS concurrent Geny processes.
# Skip-if-exists: skips samples where per-sample output already exists.
# Aggregates per-sample outputs into a single geny_hprc44.txt when done.

set -uo pipefail

PYTHON=/mypool/KIR_graph/envs/graphkir_env/bin/python3
GENY=/mypool/KIR_graph/geny/geny.py
FASTQ_ROOT=/mypool/KIR_graph/data/hprc_fastq
RESULT_DIR=/mypool/KIR_graph/results/geny
LOG_DIR=/mypool/KIR_graph/logs
AGGREGATE=$RESULT_DIR/geny_hprc44.txt
THREADS=4
PARALLEL_JOBS=10
DRY_RUN=false
SINGLE_SAMPLE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=true; shift ;;
        --sample) SINGLE_SAMPLE="$2"; shift 2 ;;
        --jobs|-j) PARALLEL_JOBS="$2"; shift 2 ;;
        --threads|-t) THREADS="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

mkdir -p "$RESULT_DIR" "$LOG_DIR"

# 44 HPRC samples: sample_id -> accession
declare -A SAMPLE_MAP=(
    [HG002]=SRR14724532
    [HG00438]=ERR3988768
    [HG005]=SRR14724518
    [HG00621]=ERR3988803
    [HG00673]=ERR3988814
    [HG00733]=ERR3988823
    [HG00735]=ERR3988824
    [HG00741]=ERR3988826
    [HG01071]=ERR3988832
    [HG01106]=ERR3988841
    [HG01109]=ERR3988842
    [HG01175]=ERR3988851
    [HG01243]=ERR3988858
    [HG01258]=ERR3988862
    [HG01358]=ERR3988875
    [HG01361]=ERR3988876
    [HG01891]=ERR3988943
    [HG01928]=ERR3988951
    [HG01952]=ERR3988958
    [HG01978]=ERR3988965
    [HG02055]=ERR3988979
    [HG02080]=ERR3988986
    [HG02109]=SRR11124438
    [HG02145]=ERR3988995
    [HG02148]=ERR3988996
    [HG02257]=ERR3989003
    [HG02572]=ERR3989028
    [HG02622]=ERR3989038
    [HG02630]=ERR3989040
    [HG02717]=ERR3989059
    [HG02723]=ERR3989060
    [HG02818]=ERR3989080
    [HG02886]=ERR3989091
    [HG03098]=ERR3989119
    [HG03453]=ERR3989166
    [HG03486]=ERR3989170
    [HG03492]=ERR3989173
    [HG03516]=ERR3989174
    [HG03540]=ERR3989177
    [HG03579]=ERR3989180
    [NA18906]=ERR3989364
    [NA19240]=ERR3989410
    [NA20129]=ERR3989456
    [NA21309]=SRR10392718
)

run_geny_sample() {
    local sample_id="$1"
    local accession="$2"
    local r1="$FASTQ_ROOT/${accession}_1.fastq.gz"
    local r2="$FASTQ_ROOT/${accession}_2.fastq.gz"
    local out="$RESULT_DIR/${sample_id}.txt"
    local log="$LOG_DIR/geny_${sample_id}.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $sample_id ($accession)" | tee -a "$log"

    if [[ -s "$out" ]]; then
        echo "  SKIP: $out already exists" | tee -a "$log"
        return 0
    fi

    if $DRY_RUN; then
        echo "  DRY-RUN: $PYTHON $GENY $r1 --file2 $r2 -t $THREADS" | tee -a "$log"
        return 0
    fi

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "  ERROR: FASTQ not found: $r1 / $r2" | tee -a "$log"
        return 1
    fi

    # Run Geny from its own directory (database.geny is a relative path in geny.py)
    local geny_dir
    geny_dir=$(dirname "$GENY")
    {
        echo "Sample: $sample_id"
        cd "$geny_dir" && "$PYTHON" "$GENY" "$r1" --file2 "$r2" -t "$THREADS" 2>>"$log"
    } > "$out"

    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        echo "  ERROR: Geny exited $exit_code" | tee -a "$log"
        rm -f "$out"
        return 1
    fi

    local lines
    lines=$(wc -l < "$out")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE $sample_id — $lines lines" | tee -a "$log"
}

export -f run_geny_sample
export PYTHON GENY FASTQ_ROOT RESULT_DIR LOG_DIR THREADS DRY_RUN

echo "=== Geny Batch Run ==="
echo "  Samples:      ${#SAMPLE_MAP[@]}"
echo "  Parallel jobs: $PARALLEL_JOBS"
echo "  Threads/job:   $THREADS"
echo "  FASTQ root:    $FASTQ_ROOT"
echo "  Result dir:    $RESULT_DIR"
echo "  Aggregate:     $AGGREGATE"
echo ""

RUNNING=0; SUCCESS=0; FAIL=0
declare -A PIDS

for sample_id in $(echo "${!SAMPLE_MAP[@]}" | tr ' ' '\n' | sort); do
    accession="${SAMPLE_MAP[$sample_id]}"

    if [[ -n "$SINGLE_SAMPLE" && "$sample_id" != "$SINGLE_SAMPLE" ]]; then
        continue
    fi

    # Wait if parallel limit reached
    while [[ $RUNNING -ge $PARALLEL_JOBS ]]; do
        wait -n 2>/dev/null || true
        for pid in "${!PIDS[@]}"; do
            if ! kill -0 "$pid" 2>/dev/null; then
                wait "$pid" 2>/dev/null
                rc=$?
                if [[ $rc -eq 0 ]]; then ((SUCCESS++)) || true
                else
                    ((FAIL++)) || true
                    echo "FAILED: ${PIDS[$pid]}" | tee -a "$LOG_DIR/geny_batch.log"
                fi
                unset "PIDS[$pid]"
                ((RUNNING--)) || true
            fi
        done
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] LAUNCH $sample_id (slot $((RUNNING+1))/$PARALLEL_JOBS)"
    run_geny_sample "$sample_id" "$accession" &
    PIDS[$!]="$sample_id"
    ((RUNNING++)) || true
done

# Drain remaining jobs
for pid in "${!PIDS[@]}"; do
    wait "$pid" 2>/dev/null
    rc=$?
    if [[ $rc -eq 0 ]]; then ((SUCCESS++)) || true
    else
        ((FAIL++)) || true
        echo "FAILED: ${PIDS[$pid]}" | tee -a "$LOG_DIR/geny_batch.log"
    fi
done

echo ""
echo "=== Batch Complete ==="
echo "  Success: $SUCCESS"
echo "  Failed:  $FAIL"

# Aggregate per-sample outputs into single file (sorted by sample_id)
echo ""
echo "=== Aggregating into $AGGREGATE ==="
> "$AGGREGATE"
for sample_id in $(echo "${!SAMPLE_MAP[@]}" | tr ' ' '\n' | sort); do
    out="$RESULT_DIR/${sample_id}.txt"
    if [[ -s "$out" ]]; then
        cat "$out" >> "$AGGREGATE"
        echo "" >> "$AGGREGATE"
    else
        echo "  WARN: missing output for $sample_id"
    fi
done

agg_lines=$(wc -l < "$AGGREGATE")
echo "  Aggregated: $agg_lines lines -> $AGGREGATE"

echo ""
echo "=== Done. To evaluate: ==="
echo "  python benchmarks/scripts/compare_geny_functional.py \\"
echo "    --geny-output $AGGREGATE"

#!/bin/bash
# ================================================================
#  Purpose: GCor pipeline
#  Author: duwenjie
# ================================================================

set -euo pipefail

# -------------------- 获取 pipeline 路径 --------------------
PIPELINE_HOME="$(cd "$(dirname "$0")" && pwd)"

# -------------------- 颜色 --------------------
RED="\033[31m"
GREEN="\033[32m"
YELLOW="\033[33m"
BLUE="\033[34m"
NC="\033[0m"

# -------------------- 时间戳日志函数 --------------------
log() {
    local level="$1"
    shift
    local msg="$*"
    local ts
    ts=$(date "+%Y-%m-%d %H:%M:%S")

    case "$level" in
        INFO)    echo -e "${BLUE}[GCor::$ts] [INFO]${NC}    $msg" ;;
        OK)      echo -e "${GREEN}[GCor::$ts] [OK]${NC}      $msg" ;;
        WARN)    echo -e "${YELLOW}[GCor::$ts] [WARN]${NC}    $msg" ;;
        ERROR)   echo -e "${RED}[GCor::$ts] [ERROR]${NC}   $msg" ;;
        *)       echo -e "[GCor::$ts] $msg" ;;
    esac

    # 写入日志文件（如果已定义）
    [[ -n "${LOG_FILE:-}" ]] && echo "[GCor::$ts] [$level] $msg" >> "$LOG_FILE"
}

# -------------------- 参数帮助 --------------------
usage() {
cat <<EOF
Usage: $(basename "$0") [options]

Run Hi-C full pipeline:
  1. Convert BAM to .mcool
  2. Process matrix
  3. Detect abnormal trend
  4. Analyse abnormal trend
  5. Locate breakpoint

Required arguments:
  -f, --fasta     <fasta>          Reference FASTA file
  -b, --bam       <bam/pairs>      Input Hi-C BAM file
  -o, --output    <prefix>         Output prefix (also used for log file)

Optional arguments:
  --run-trf                        Run TRF (Tandem Repeat Finder) before matrix processing
  --diag-band-width <bins>         Half-width of diagonal bands (default: 0)
  --trend-p <float>                Row/column trend p-value threshold (default: 0.05)
  --structure-p <float>            Structural feature p-value threshold (default: 0.05)
  --structure-z <float>            Structural feature absolute Z threshold (default: 1.96)

Example:
  $(basename "$0") -f ref.fa -b hic.bam -o sample
EOF
exit 0
}

# -------------------- 参数解析 --------------------
FA_FILE=""
BAM_FILE=""
OUT_PREFIX=""
RUN_TRF=false  # 默认不运行 TRF
DIAG_BAND_WIDTH=0
TREND_P=0.05
STRUCTURE_P=0.05
STRUCTURE_Z=1.96

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--fasta)   FA_FILE="$2"; shift 2 ;;
        -b|--bam)     BAM_FILE="$2"; shift 2 ;;
        -o|--output)  OUT_PREFIX="$2"; shift 2 ;;
        --run-trf)    RUN_TRF=true; shift 1 ;; 
        --diag-band-width) DIAG_BAND_WIDTH="$2"; shift 2 ;;
        --trend-p) TREND_P="$2"; shift 2 ;;
        --structure-p) STRUCTURE_P="$2"; shift 2 ;;
        --structure-z) STRUCTURE_Z="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *)            log ERROR "Unknown option: $1"; usage ;;
    esac
done

# -------------------- 检查参数 --------------------
[[ -z "$FA_FILE" || -z "$BAM_FILE" || -z "$OUT_PREFIX" ]] && {
    log ERROR "Missing required parameters."
    usage
}

[[ ! -f "$FA_FILE" ]] && { log ERROR "FASTA not found: $FA_FILE"; exit 1; }
[[ ! -f "$BAM_FILE" ]] && { log ERROR "BAM not found: $BAM_FILE"; exit 1; }

# -------------------- 初始化日志文件 --------------------
LOG_FILE="${OUT_PREFIX}.log"
echo "================ Pipeline Log ================" > "$LOG_FILE"
echo "START: $(date)" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

log INFO "Pipeline started."
log INFO "FASTA:  $FA_FILE"
log INFO "BAM:    $BAM_FILE"
log INFO "Prefix: $OUT_PREFIX"
log OK   "Logging to: $LOG_FILE"


# -------------------- Step 1: BAM → mcool --------------------
log INFO "Step 1: Convert BAM to mcool"
if bash "${PIPELINE_HOME}/bam2mcool.sh" -i "$BAM_FILE" -r "$FA_FILE" -p "$OUT_PREFIX" >>"$LOG_FILE" 2>&1; then
    log OK "mcool created successfully."
else
    log ERROR "bam2mcool.sh failed. See $LOG_FILE"
    exit 1
fi

# -------------------- Step 2: Matrix processing --------------------
if [[ "$RUN_TRF" == true ]]; then
    log INFO "Step 2: Process matrix"
    if bash "${PIPELINE_HOME}/process.sh" -f "$FA_FILE" -m "${OUT_PREFIX}.mcool" -o "$OUT_PREFIX" -T >>"$LOG_FILE" 2>&1; then
        log OK "Matrix processing complete."
    else
        log ERROR "process.sh failed."
        exit 1
    fi
else
    log INFO "Step 2: Process matrix(TRF skipped)"
    if bash "${PIPELINE_HOME}/process.sh" -f "$FA_FILE" -m "${OUT_PREFIX}.mcool" -o "$OUT_PREFIX" >>"$LOG_FILE" 2>&1; then
        log OK "Matrix processing complete."
    else
        log ERROR "process.sh failed."
        exit 1
    fi
fi

# -------------------- Step 3: Abnormal trend detection --------------------
if [[ "$RUN_TRF" == true ]]; then
    log INFO "Step 3: Detect abnormal trends"
    if python "${PIPELINE_HOME}/find_abnormal_trend_from_matrix.py" \
        --mcool_file "${OUT_PREFIX}.mcool" \
        --boundary_file "${OUT_PREFIX}_100000_boundaries.gff" \
        --tr_file "${OUT_PREFIX}_TR_F2000_M100000.txt" \
        --gap_file "${OUT_PREFIX}_GAP.txt" \
        --diag_band_width "$DIAG_BAND_WIDTH" >>"$LOG_FILE" 2>&1; then
        log OK "Abnormal trend detection complete."
    else
        log ERROR "find_abnormal_trend_from_matrix.py failed."
        exit 1
    fi
else
    log INFO "Step 3: Detect abnormal trends"
    if python "${PIPELINE_HOME}/find_abnormal_trend_from_matrix.py" \
        --mcool_file "${OUT_PREFIX}.mcool" \
        --boundary_file "${OUT_PREFIX}_100000_boundaries.gff" \
        --gap_file "${OUT_PREFIX}_GAP.txt" \
        --diag_band_width "$DIAG_BAND_WIDTH" >>"$LOG_FILE" 2>&1; then
        log OK "Abnormal trend detection complete."
    else
        log ERROR "find_abnormal_trend_from_matrix.py failed."
        exit 1
    fi
fi

# -------------------- Step 4: Per-chromosome analysis --------------------
log INFO "Step 4: Per-chromosome submatrix analysis"

RESULT_DIR="results"
if [[ ! -d "$RESULT_DIR" ]]; then
    log ERROR "Directory 'results/' not found!"
    exit 1
fi

for chr_dir in ${RESULT_DIR}/*/; do
    log INFO "Processing: $chr_dir"

    sub="${chr_dir}/submatrix_trend_data"
    if [[ ! -d "$sub" ]]; then
        log WARN "Missing submatrix_trend_data/, skipped."
        continue
    fi

    pushd "$sub" >/dev/null

    python "${PIPELINE_HOME}/analyse_abnormal_trend.py" --input *.tsv \
        --trend_p "$TREND_P" \
        --structure_p "$STRUCTURE_P" \
        --structure_z "$STRUCTURE_Z" >>"$LOG_FILE" 2>&1
    python "${PIPELINE_HOME}/locate_breakpoint_position.py" \
        -i "$(realpath ../../../$(basename "$BAM_FILE"))" \
        -r result_cluster_merged_coords.tsv >>"$LOG_FILE" 2>&1

    popd >/dev/null
done

log OK "Pipeline finished successfully!"
echo "END: $(date)" >> "$LOG_FILE"

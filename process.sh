#!/bin/bash
# ================================================================
#  Script: process.sh
#  Author: <duwenjie>
#  Purpose: Run TRF, identify gaps, and call TADs from a .mcool file
#  Dependencies: trf, seqkit, hicexplorer
# ================================================================

set -euo pipefail


# -------------------- 参数帮助 --------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Run TRF, locate gaps, and call TADs from Hi-C .mcool data.

Required arguments:
  -f <fasta>           Input genome FASTA file
  -m <mcool>           Input Hi-C multi-resolution .mcool file
  -o <prefix>          Output prefix for all results

Optional arguments:
  -r <resolution>      Resolution to use for TAD calling (default: 100000)
  -t <threads>         Number of threads for hicFindTADs (default: 12)
  -h                   Show this help message and exit

TRF result arguments:
  -T                   ENABLE TRF (default: disabled)
  -P <TRF_N_PROCS>     Parallel jobs (default: 4)
  -F <tr_filter_len>   Minimum TR length to keep (default: 2000)
  -M <tr_merge_len>    Minimum merge distance between TRs (default: 100000)

Example:
  $(basename "$0") -f genome.fa -m sample.mcool -o sample -r 100000 -t 12

Output:
  sample_TRF.txt       # Tandem repeats output
  sample_GAP.txt       # Regions with long Ns
  sample_<RES>_domains.txt # TAD calls at specified resolution
EOF
    exit 0
}

# -------------------- 参数解析 --------------------
FASTA=""
MCOOL=""
OUTPUT_PREFIX=""
RES="100000"
THREADS=12

ENABLE_TRF=0 
TRF_N_PROCS=4
TR_FILTER_LEN=2000
TR_MERGE_LEN=100000

while getopts "f:m:o:r:t:h" opt; do
    case $opt in
        f) FASTA=$OPTARG ;;
        m) MCOOL=$OPTARG ;;
        o) OUTPUT_PREFIX=$OPTARG ;;
        r) RES=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        T) ENABLE_TRF=1 ;; 
        P) TRF_N_PROCS=$OPTARG ;;
        F) TR_FILTER_LEN=$OPTARG ;;
        M) TR_MERGE_LEN=$OPTARG ;;
        h) usage ;;
        *) usage ;;
    esac
done

# -------------------- 参数检查 --------------------
if [[ -z "$FASTA" || -z "$MCOOL" || -z "$OUTPUT_PREFIX" || -z "$RES" ]]; then
    echo "[process::Error] Missing required argument(s)."
    usage
fi

# -------------------- 日志函数 --------------------
log() {
    echo "[process::$(date '+%F %T')] $*"
}

log "Start processing"
log "FASTA: $FASTA"
log "mcool: $MCOOL"
log "Output prefix: $OUTPUT_PREFIX"
log "Resolution: $RES"
log "Threads: $THREADS"
log "TRF running processes: $TRF_N_PROCS"
log "TR filter length: $TR_FILTER_LEN"
log "TR merge distance: $TR_MERGE_LEN"
echo

# ================================================================
#   Step 1: TRF (only if ENABLE_TRF=1)
# ================================================================
if [[ "$ENABLE_TRF" -eq 1 ]]; then
    log "TRF enabled. Running tandem repeat annotation."

    TMP_DIR=$(mktemp -d)
    trap 'rm -rf "$TMP_DIR"' EXIT  
    seqkit split --by-id -O "$TMP_DIR" "$FASTA"

    FINAL_TR_RAW="${OUTPUT_PREFIX}_TR.txt"
    FINAL_TR_FILT="${OUTPUT_PREFIX}_TR_F${TR_FILTER_LEN}_M${TR_MERGE_LEN}.txt"
    : > "$FINAL_TR_RAW"
    : > "$FINAL_TR_FILT"

    process_chr() {
        local chr_fa="$1"
        local base=$(basename "$chr_fa")
        chr_tmp=${base%.fa}
        chr_name="${chr_tmp##*_}"

        trf "$chr_fa" 2 7 7 80 10 50 500 -d -h 
        local dat="${base}.2.7.7.80.10.50.500.dat"
        if [[ ! -f "$dat" ]]; then
            echo "[CHR:$chr_name] TRF failed!" >&2
            return 1
        fi

        sed -i "1i #CHR:${chr_name}" "$dat"

        awk -v chr="$chr_name" \
            -v flen="$TR_FILTER_LEN" \
            -v mlen="$TR_MERGE_LEN" \
            -v raw="$FINAL_TR_RAW" \
            -v filt="$FINAL_TR_FILT" \
            -v last=0 \
            'BEGIN {OFS="\t"}
            /^#/ {print > raw; next}
            {
                if ($1 ~ /^[0-9]+$/) print chr, $0 > raw;
                if (($2-$1)>=flen && ($1-last)>=mlen) {
                    print chr, $0 > filt; last=$1;
                }
            }' "$dat"
    }

    export -f process_chr
    PARALLEL=$(which parallel)
    find "$TMP_DIR" -name "*.fa" -print0 |
        $PARALLEL -0 -j "$TRF_N_PROCS" --no-notice process_chr {}

    log "TRF finished:"
    log "  Raw TR: $FINAL_TR_RAW"
    log "  Filtered TR: $FINAL_TR_FILT"
else
    log "TRF disabled. Skipping TRF step."
fi


# -------------------- Step 2: Identify gaps --------------------
log "Locating long N regions (gaps)"
GAP_FILE="${OUTPUT_PREFIX}_GAP.txt"
seqkit locate -r -p '"N{100,}"' "$FASTA" > "${GAP_FILE}"
awk '$4 ~ /^\+/{print $1"\t"$5}' "${GAP_FILE}" > "${GAP_FILE}.tmp" && mv "${GAP_FILE}.tmp" "${GAP_FILE}"
log "Gaps output saved to "${GAP_FILE}""

# -------------------- Step 3: Call TADs --------------------
log "Running hicFindTADs at resolution $RES"
log "Detecting chromosomes from $MCOOL ..."

hicFindTADs \
    -m "${MCOOL}::/resolutions/${RES}" \
    --outPrefix "${OUTPUT_PREFIX}_${RES}" \
    --numberOfProcessors "$THREADS" \
    --minBoundaryDistance "$RES" \
    --correctForMultipleTesting fdr \
    --thresholdComparisons 0.05

# for chr in $(cut -f1 ${FASTA}); do
#     log "Running hicFindTADs in ${chr}"
#     hicFindTADs \
#         -m "${MCOOL}::/resolutions/${RES}/${chr}" \
#         --outPrefix "${OUTPUT_PREFIX}_${chr}_${RES}" \
#         --numberOfProcessors "${THREADS}" \
#         --minBoundaryDistance "${RES}" \
#         --correctForMultipleTesting fdr \
#         --thresholdComparisons 0.05
# done


log "TAD calling finished. Output prefix: ${OUTPUT_PREFIX}_${RES}"
log "All done!"

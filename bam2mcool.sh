#!/bin/bash
# ================================================================
#  Script: bam_to_mcool.sh
#  Author: duwenjie
#  Purpose: Convert Hi-C BAM or PAIRS file into multi-resolution .mcool
#  Dependencies: samtools, pairtools, cooler
# ================================================================

set -euo pipefail

# -------------------- 参数帮助 --------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Convert Hi-C BAM or PAIRS file to multi-resolution .mcool matrix.

Required arguments:
  -i <input>        Input file (.bam or .pairs/.pairs.gz)
  -r <reference>    Reference genome FASTA (used to generate chrom.sizes)
  -p <prefix>       Output prefix name (used for all intermediate/output files)

Optional arguments:
  -b <binsize>      Base bin size (default: 10000)
  -q <mapq>         Minimum MAPQ filter threshold (default: 0)
  -t <threads>      Number of parallel threads (default: 8)
  -B <yes|no>       Whether to perform matrix balancing (default: yes)
  -h                Show this help message and exit

Example:
  $(basename "$0") -i aligned.bam -r ref.fa -p sample -b 25000 -q 20 -t 8 -B yes

Output:
  sample.mcool      # multi-resolution contact matrix
  sample_10000.cool # base-resolution matrix before zoomify
  sample.chrom.sizes
EOF
    exit 0
}

# -------------------- 参数解析 --------------------
INPUT=""
REF=""
PREFIX=""
BINSIZE=""
MAPQ=""
THREADS=""
DO_BALANCE=""

while getopts "i:r:p:b:q:t:B:h" opt; do
  case $opt in
    i) INPUT=$OPTARG ;;
    r) REF=$OPTARG ;;
    p) PREFIX=$OPTARG ;;
    b) BINSIZE=$OPTARG ;;
    q) MAPQ=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    B) DO_BALANCE=$OPTARG ;;   
    h) usage ;;
    *) usage ;;
  esac
done

# -------------------- 参数检查 --------------------
if [[ -z "$INPUT" || -z "$REF" || -z "$PREFIX" ]]; then
  echo "[bam2mcool::Error] Missing required argument(s)."
  usage
fi

# -------------------- 默认参数（仅当未输入时） --------------------
BINSIZE=${BINSIZE:-10000}
MAPQ=${MAPQ:-0}
THREADS=${THREADS:-21}
DO_BALANCE=${DO_BALANCE:-yes}   # ✅ 默认执行平衡

# -------------------- 输出文件路径 --------------------
CHROMSIZES="${PREFIX}.chrom.sizes"
PAIRS_RAW="${PREFIX}.pairs"
COOL="${PREFIX}_${BINSIZE}.cool"
MCOOL="${PREFIX}.mcool"

# -------------------- 日志函数 --------------------
log() {
  echo "[bam2mcool::`date '+%F %T'`] $*"
}

# -------------------- 开始流程 --------------------
log "Start Hi-C matrix generation pipeline"
log "Input: $INPUT"
log "Reference: $REF"
log "Prefix: $PREFIX"
log "Bin size: $BINSIZE bp"
log "MAPQ threshold: $MAPQ"
log "Threads: $THREADS"
log "Do balance: $DO_BALANCE"

## Step 1. 生成染色体长度文件
if [ ! -f "$CHROMSIZES" ]; then
  log "Generating chrom.sizes from $REF"
  samtools faidx "$REF"
  cut -f1,2 "${REF}.fai" > "$CHROMSIZES"
fi

## Step 2. 判断输入类型
EXT="${INPUT##*.}"
if [[ "$EXT" == "bam" ]]; then
  log "Input detected as BAM"
  log "Parsing BAM -> .pairs.gz"
  pairtools parse -c "$CHROMSIZES" \
    --add-columns mapq \
    --nproc-in $THREADS --nproc-out $THREADS \
    -o "$PAIRS_RAW" "$INPUT"
elif [[ "$EXT" == "pairs" || "$EXT" == "gz" ]]; then
  log "Input detected as PAIRS"
  cp "$INPUT" "$PAIRS_RAW"
else
  echo "[Error] Unsupported input format: $INPUT"
  exit 1
fi

## Step 3. check header
PAIRS_RAW_FILE="$PAIRS_RAW"
if [[ -f "$PAIRS_RAW_FILE" ]]; then
    if ! grep -q "^#chromsize:" "$PAIRS_RAW_FILE"; then
        log "No header found in $PAIRS_RAW_FILE, adding standard header from chrom.sizes"
        tmp="${PAIRS_RAW_FILE%.pairs}_with_header.pairs"
        {
            echo "## pairs format v1.0.0"
            echo "#shape: upper triangle"
            echo "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"
            awk '{printf("#chromsize: %s %s\n",$1,$2)}' "$CHROMSIZES"
            cat "$PAIRS_RAW_FILE"
        }  > "$tmp"
        PAIRS_RAW="$tmp"
    fi
else
    echo "[Error] PAIRS file $PAIRS_RAW_FILE not found!"
    exit 1
fi


## Step 4. 压缩pairs
log "Compress pairs"
bgzip -c "$PAIRS_RAW" > "$PAIRS_RAW".gz

## Step 5. 排序pairs.gz
log "Sort pairs.gz"
pairtools sort --nproc $THREADS -o ${PAIRS_RAW%.pairs}.sorted.pairs.gz "$PAIRS_RAW".gz

# Step 6. 构建索引
log "Build index sorted.pairs.gz"
pairix -s 2 -b 3 -d 4 -u 5 ${PAIRS_RAW%.pairs}.sorted.pairs.gz

# Step 7. 生成单分辨率 .cool
log "Creating ${BINSIZE}-bin .cool"
cooler cload pairix "${CHROMSIZES}:${BINSIZE}" ${PAIRS_RAW%.pairs}.sorted.pairs.gz "$COOL"

# Step 8. 平衡矩阵（可选）
if [[ "$DO_BALANCE" == "yes" ]]; then
  log "Balancing .cool matrix"
  cooler balance --nproc $THREADS "$COOL"
else
  log "Skipping matrix balancing as requested"
fi

# Step 9. 生成多分辨率 .mcool
log "Zoomifying to .mcool"
if [[ "$DO_BALANCE" == "yes" ]]; then
  cooler zoomify --balance --nproc $THREADS "$COOL" -o "$MCOOL" -r N
else
  cooler zoomify --nproc $THREADS "$COOL" -o "$MCOOL" -r N
fi

if [[ -n "$tmp" ]]; then
    rm "$tmp"
fi
log "All done!"

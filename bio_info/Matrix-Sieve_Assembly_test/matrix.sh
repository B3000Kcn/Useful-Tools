# ... (前面 seqkit fx2tab 部分保持不变) ...

echo "Calculating Adaptive Thresholds (Pure Bash/Awk)..."

# 1. 计算总行数 (Total Reads)
TOTAL_READS=$(wc -l < read_metrics.tsv)
echo "Total Reads: ${TOTAL_READS}"

# 2. 定义计算函数：输入列号(2=Len, 3=Qual)，输入百分比(0.05=95%留存)
# 逻辑：排序 -> 计算目标行号 -> 提取该行数值
get_quantile() {
    col=$1
    p=$2
    # 计算目标行号 (Total * P)，取整
    target_line=$(awk -v t="${TOTAL_READS}" -v p="${p}" 'BEGIN {printf "%.0f", t * p}')
    # 如果算出来是0，强制设为1
    [[ "$target_line" -eq 0 ]] && target_line=1
    
    # 核心：对指定列排序，取第 target_line 行的数值
    # sort -k${col},${col}n : 按指定列数值排序 (从小到大)
    # head -n : 取前n行
    # tail -n 1 : 取最后一行(即第n行)
    # cut -f : 提取该列
    sort -t$'\t' -k${col},${col}n read_metrics.tsv | head -n "${target_line}" | tail -n 1 | cut -f${col}
}

# 3. 计算 Length 阈值 (第2列)
# 5%分位数 = 丢弃最烂的5% = 留存95%
L95=$(get_quantile 2 0.05)
L80=$(get_quantile 2 0.20)
L60=$(get_quantile 2 0.40)

# 4. 计算 Quality 阈值 (第3列)
Q95=$(get_quantile 3 0.05)
Q80=$(get_quantile 3 0.20)
Q60=$(get_quantile 3 0.40)

echo "Adaptive L-Grids (95%/80%/60%): ${L95}, ${L80}, ${L60}"
echo "Adaptive Q-Grids (95%/80%/60%): ${Q95}, ${Q80}, ${Q60}"

# ... (后面生成任务部分保持不变) ...

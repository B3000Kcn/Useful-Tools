## 🧬 丝滑小连招：完美单序列抽离术 (Single Sequence Extraction)

**场景**：从巨大的全基因组比对结果中，干净地剥离出某一条特定序列（如质粒、线粒体）的比对数据，并彻底清洗 Header，防止在 IGV/Qualimap 中出现数百条空轨道。

### 🕹️ 连招一：比对与生成 SAM (Mapping)
*注：显式指定绝对路径，解决 minimap2 与 samtools 环境冲突问题。*

```bash
# 建立索引并比对
/path/to/minimap2 -d /path/to/ref.mmi /path/to/ref.fa && \
/path/to/minimap2 -t 128 -ax map-ont /path/to/ref.fa /path/to/reads.fq.gz > /path/to/output.sam
```
> **解析**：
> *   `/path/to/minimap2`：填入你的 minimap2 绝对路径。
> *   `-ax map-ont`：针对 ONT 数据的比对参数。

### 🕹️ 连招二：转 BAM、排序与清理 (Sort & Clean)
*注：转换后立即删除巨大的 SAM 文件，释放空间。*

```bash
/path/to/samtools view -@ 128 -Sb /path/to/output.sam > /path/to/output.bam && \
/path/to/samtools sort -@ 12 /path/to/output.bam -o /path/to/output_sorted.bam && \
rm -rf /path/to/output.bam /path/to/output.sam
```
> **解析**：
> *   `-Sb`：将 SAM 转为 BAM。
> *   `rm -rf ...`：**关键动作**，生成的中间文件巨大，用完即焚。

### 🕹️ 连招三：究极抽离与 Header 清洗 (Extraction & Header Cleaning)
*注：这是整套连招的灵魂。`grep` 负责“外科手术式”地切除无关序列的 Header 定义。*

```bash
# 请将下方的 "Target_Seq_ID" 替换为你实际要抽离的序列名（如 "2" 或 "plasmid_A"）
/path/to/samtools index /path/to/output_sorted.bam && \
/path/to/samtools view -h /path/to/output_sorted.bam "Target_Seq_ID" | \
grep -E "^@HD|^@PG|^@CO|SN:Target_Seq_ID|^[^@]" | \
/path/to/samtools view -b - > /path/to/final_extracted.bam
```
> **解析**：
> *   `view ... "Target_Seq_ID"`：只提取比对到该序列的 Reads。
> *   `grep -E ...`：
>     *   `^@HD|^@PG|^@CO`：保留通用头信息。
>     *   **`SN:Target_Seq_ID`**：**核心！只保留目标序列的定义行 (@SQ)，丢弃其他成百上千条无关序列的定义。**
>     *   `^[^@]`：保留正文（比对记录）。

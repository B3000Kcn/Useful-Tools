# BLAST DB Parallel Extract FASTA

利用 `GNU parallel` 对 NCBI 多分卷 BLAST 数据库进行高并发、流式无中间文件全量提取 FASTA 的高性能 POSIX Shell 工具。

---

## 核心特性

- **纯内存流式合并**：通过字符序列化将每条 FASTA 记录（序列名与序列正文）压扁为传输单行，借助管道直接汇入最终文件，全程无需在磁盘生成任何分卷中间文件，零额外写放大。
- **原子级防混写与防错位**：序列名与序列正文在子进程内即通过 `\0`（NUL byte）紧密绑定为单一传输行，配合 `parallel --line-buffer` 的行级原子流转，彻底杜绝并发写入时的数据混写、断行或序列名错位。
- **Fail-Fast 异常熔断**：管道末端强制校验每条记录的引导符（`die unless /^>/`），并在捕获上游异常标记时即刻中止，确保绝不静默生成残缺或损坏数据。
- **动态别名补齐与排他保护**：自动扫描底层的 `.nsq` 或 `.psq` 物理数据卷并按需建立别名（退出自动清理）；内置目录级原子排他锁，防止多任务并发冲突。
- **纯血 POSIX 规范**：以 `#!/bin/sh` 与 `set -eu` 编写，无 Bashism 专有依赖。

---

## 依赖环境

使用前请确保以下工具位于系统环境变量 `$PATH` 中：

- `blastdbcmd` (来自 NCBI BLAST+ 软件包)
- `parallel` (GNU parallel)
- `perl` (v5.10 及以上)
- 基础系统工具：`sh`, `find`, `sed`, `sort`, `mkdir`, `rm`

---

## 参数说明

```text
用法: sh blast_db_parallel_extract_fasta.sh -n|-p --db-dir <DIR> --db-prefix <PREFIX> --output <PATH> [选项]

必需参数:
  -n, --nucl            核酸数据库 (Nucleotide)
  -p, --prot            蛋白数据库 (Protein)
                        (注: 必须显式指定 -n 或 -p 之一)
  --db-dir DIR          BLAST 数据库所在目录
  --db-prefix PREFIX    BLAST 数据库公共前缀 (例如: nt, nr)
  --output PATH         输出目标路径:
                        - 默认模式: 最终合并生成的 FASTA 文件路径 (如 /path/to/all.fasta)
                        - 分卷模式: 存放各分卷 .fasta 文件的目标目录路径

可选参数:
  -j, --jobs INT        并行任务数 (默认: 自动检测宿主机可用逻辑核心数)
  --split               分卷输出模式: 为每个分卷生成独立的 .fasta 文件
  -h, --help            显示帮助信息并退出
```

---

## 使用示例

### 1. 全量流式合并（默认推荐）

提取多分卷核酸数据库（如 `nt`）并直接合并为单一 FASTA 文件，零中间文件：

```bash
sh blast_db_parallel_extract_fasta.sh \
    -n \
    --db-dir /path/to/blast_db \
    --db-prefix nt \
    --output /path/to/output/nt.all.fasta \
    --jobs 64
```

提取多分卷蛋白数据库（如 `nr`）：

```bash
sh blast_db_parallel_extract_fasta.sh \
    -p \
    --db-dir /path/to/blast_db \
    --db-prefix nr \
    --output /path/to/output/nr.all.fasta \
    --jobs 32
```

### 2. 分卷导出模式 (`--split`)

若需要单独保留各分卷独立的 `.fasta` 文件：

```bash
sh blast_db_parallel_extract_fasta.sh \
    -n \
    --db-dir /path/to/blast_db \
    --db-prefix nt \
    --output /path/to/output_dir \
    --jobs 16 \
    --split
```

---

## 工作机制

```text
[分卷 00] -> blastdbcmd -> Perl: 序列内部与名后 \n 替换为 \0 ---\
[分卷 01] -> blastdbcmd -> Perl: 序列内部与名后 \n 替换为 \0 ----+--> GNU parallel --line-buffer
   ...                                                           |    (整行原子无锁汇入主管道)
[分卷 N ] -> blastdbcmd -> Perl: 序列内部与名后 \n 替换为 \0 ---/                 |
                                                                                    v
                                                                        Perl: tr/\0/\n/ 查表级内存还原
                                                                              die unless /^>/;
                                                                                    |
                                                                                    v
                                                                           输出最终合并 FASTA
```

1. **序列名与序列原子化绑定**：Worker 端 Perl 在读取 `blastdbcmd` 时，将序列名之后及序列内部的所有换行统一置换为 `\0`（NUL 字节），使**整条 FASTA 记录在传输中成为唯一的逻辑单行**；
2. **无锁原子流转**：通过 `parallel --line-buffer` 针对单行记录的原子刷新特性，各进程数据完整入流，不会相互插入被打断；
3. **查表还原与熔断**：主管道末端以 Perl `tr/\0/\n/` 执行纳秒级内存查表还原；若任意数据行未以 `>` 开头或捕获到错误标记，即刻触发 `die` 熔断中断管道。

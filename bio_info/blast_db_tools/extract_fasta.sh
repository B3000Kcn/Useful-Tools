#!/bin/sh
set -eu

usage() {
    cat <<EOF
用法：
  合并：sh $0 --db-dir DIR --db-prefix PREFIX --output FILE [--jobs N]
  分卷：sh $0 --db-dir DIR --db-prefix PREFIX --output DIR  [--jobs N] --split

必填参数：
  --db-dir DIR         BLAST 数据库目录
  --db-prefix PREFIX  数据库公共前缀，例如 nt_prok
  --output PATH       合并模式：大 FASTA 的完整路径；分卷模式：输出目录

可选参数：
  --jobs N, -j N      并行进程数，正整数，默认 32
  --split             输出分卷 FASTA；不加此参数则直接合并
  --help, -h          显示帮助

注意：同名输出文件会被覆盖；合并结果不保证分卷顺序。
EOF
}

die() { printf '错误：%s\n' "$*" >&2; exit 1; }
need_value() {
    [ "$#" -ge 2 ] && [ -n "$2" ] || die "$1 后面必须提供参数值。"
    case "$2" in --*|-j|-h) die "$1 后面必须提供参数值。" ;; esac
}

# --- 命令行参数；只有并行数和输出模式设有默认值 ---
DB_DIR=
DB_PREFIX=
OUTPUT=
NUM_JOBS=32
MERGE_FASTA=1

while [ "$#" -gt 0 ]; do
    case "$1" in
        --db-dir)    need_value "$@"; DB_DIR=$2; shift 2 ;;
        --db-prefix) need_value "$@"; DB_PREFIX=$2; shift 2 ;;
        --output)    need_value "$@"; OUTPUT=$2; shift 2 ;;
        --jobs|-j)   need_value "$@"; NUM_JOBS=$2; shift 2 ;;
        --split)     MERGE_FASTA=0; shift ;;
        --help|-h)   usage; exit 0 ;;
        *)           die "未知参数：$1。使用 --help 查看用法。" ;;
    esac
done

# 缺少任何必填参数均退出，不启动提取、不创建输出。
[ -n "$DB_DIR" ] || die "缺少必填参数：--db-dir"
[ -n "$DB_PREFIX" ] || die "缺少必填参数：--db-prefix"
[ -n "$OUTPUT" ] || die "缺少必填参数：--output"
case "$NUM_JOBS" in
    ''|*[!0-9]*|0*) die "--jobs 必须是正整数，例如 32。" ;;
esac
case "$DB_PREFIX" in
    *[!a-zA-Z0-9_.-]*) die "--db-prefix 只允许字母、数字、下划线、点和连字符。" ;;
esac

# 将相对路径转换为绝对路径，便于子进程使用。
case "$DB_DIR" in /*) ;; *) DB_DIR="$PWD/$DB_DIR" ;; esac
case "$OUTPUT" in /*) ;; *) OUTPUT="$PWD/$OUTPUT" ;; esac
[ -d "$DB_DIR" ] || die "数据库目录不存在：$DB_DIR"
DB_DIR=$(CDPATH= cd "$DB_DIR" && pwd -P) || die "无法进入数据库目录。"

for cmd in blastdbcmd parallel perl find sed sort dirname mkdir; do
    command -v "$cmd" >/dev/null 2>&1 || die "找不到命令：$cmd"
done

# 提前检查是否存在匹配分卷，避免没有输入时创建空结果。
found=0
for nal in "$DB_DIR"/"$DB_PREFIX".*.nal; do
    [ -f "$nal" ] || continue
    found=1
    break
done
[ "$found" = 1 ] || die "未找到 $DB_DIR/${DB_PREFIX}.*.nal"

if [ "$MERGE_FASTA" = 1 ]; then
    [ ! -d "$OUTPUT" ] || die "合并模式的 --output 必须是文件路径，不能是目录。"
    FASTA_OUTPUT_DIR=$(dirname "$OUTPUT")
    OUTPUT_TARGET=$OUTPUT
else
    FASTA_OUTPUT_DIR=$OUTPUT
    OUTPUT_TARGET=/dev/null
fi

# 每条 FASTA 记录封装成一个传输行，避免并行时标题与序列混写。
# 仅经内存管道传输；末端恢复换行，不创建分卷中间文件。
PACK_FASTA='
    open(my $in, "-|", "blastdbcmd", "-db", $ARGV[0], "-entry", "all")
        or die "Cannot run blastdbcmd: $!\n";
    my $started = 0;
    while (my $line = <$in>) {
        chomp $line;
        if ($line =~ /^>/) {
            print "\n" if $started;
            $started = 1;
        } else {
            print "\0";
        }
        print $line;
    }
    print "\n" if $started;
    close($in) or exit 1;
'
export DB_DIR FASTA_OUTPUT_DIR MERGE_FASTA PACK_FASTA
PARALLEL_SHELL=/bin/sh
export PARALLEL_SHELL

# --- 脚本主逻辑 ---
echo "并行提取 FASTA 脚本开始执行；并行进程数：$NUM_JOBS" >&2

if [ ! -d "$FASTA_OUTPUT_DIR" ]; then
    echo "正在创建输出目录：$FASTA_OUTPUT_DIR" >&2
    mkdir -p "$FASTA_OUTPUT_DIR" || exit 1
fi

if [ "$MERGE_FASTA" = 1 ]; then
    echo "直接合并到：$OUTPUT" >&2
else
    echo "分卷输出目录：$OUTPUT" >&2
fi

# 保留 find -> sed -> sort -> parallel 的处理结构。
find "$DB_DIR" -maxdepth 1 -type f -name "${DB_PREFIX}.*.nal" \
    | sed 's|.*/||; s/\.nal$//' \
    | sort \
    | {
        parallel --plain -j "$NUM_JOBS" --eta --line-buffer --halt soon,fail=1 '
            db_name={}
            db_path="$DB_DIR/$db_name"
            output_path="$FASTA_OUTPUT_DIR/${db_name}.fasta"

            echo "正在处理分卷：$db_name ..." >&2

            if [ "$MERGE_FASTA" = 1 ]; then
                perl -e "$PACK_FASTA" "$db_path" || exit 1
            else
                blastdbcmd -db "$db_path" -entry all -out "$output_path" || exit 1
            fi

            echo "分卷 $db_name 处理完成。" >&2
        ' || printf '\n__EXTRACTION_FAILED__\n'
    } \
    | perl -ne '
        die "提取失败，请检查上方日志。\n" unless /^>/;
        tr/\0/\n/;
        print or die "写入失败：$!\n";
        END {
            die "没有提取到序列。\n" if $ENV{MERGE_FASTA} && !$.;
            close STDOUT or die "写入失败：$!\n";
        }
    ' > "$OUTPUT_TARGET" || {
        echo "执行失败，已写入的输出可能不完整。" >&2
        exit 1
    }

echo "所有分卷均已处理完毕！输出：$OUTPUT" >&2

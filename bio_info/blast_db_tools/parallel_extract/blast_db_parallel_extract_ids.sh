#!/bin/sh
set -eu

usage() {
    cat <<EOF
用法：
  合并：sh $0 (-n | -p) --db-dir DIR --db-prefix PREFIX --output FILE [--jobs N]
  分卷：sh $0 (-n | -p) --db-dir DIR --db-prefix PREFIX --output DIR  [--jobs N] --split

必填参数：
  -n 或 -p            必须二选一：-n 核酸库，-p 蛋白库；没有默认类型
  --db-dir DIR         BLAST 数据库目录
  --db-prefix PREFIX  数据库公共前缀，例如 nt_prok
  --output PATH       合并模式：总 ID 文件路径；分卷模式：输出目录

可选参数：
  --jobs N, -j N      并行进程数，正整数，默认 32
  --split             输出分卷 ID 文件；不加此参数则直接合并
  --help, -h          显示帮助

注意：沿用原脚本的 -outfmt "%a"；同名输出文件会被覆盖。
合并结果不保证分卷顺序，不去重、不排序 ID。
别名处理：-n 扫描 PREFIX.*.nsq 并补 .nal；-p 扫描 PREFIX.*.psq 并补 .pal。
已有非空别名不改写；缺失或空的别名自动补齐。
退出时删除本次分卷的 .nal/.pal（包括原已存在的）；保留顶层总库别名。
同一数据库目录、同一前缀的两份脚本请依次运行。
不能清理 kill -9、断电等无法捕获的强制中断留下的文件。

EOF
}

die() { printf '错误：%s\n' "$*" >&2; exit 1; }
need_value() {
    [ "$#" -ge 2 ] && [ -n "$2" ] || die "$1 后面必须提供参数值。"
    case "$2" in -*) die "$1 后面必须提供参数值。" ;; esac
}

# --- 命令行参数：仅并行数和输出模式设有默认值 ---
DB_TYPE=
DB_DIR=
DB_PREFIX=
OUTPUT=
NUM_JOBS=32
MERGE_IDS=1

while [ "$#" -gt 0 ]; do
    case "$1" in
        -n) [ "$DB_TYPE" != prot ] || die "-n 和 -p 不能同时指定。"; DB_TYPE=nucl; shift ;;
        -p) [ "$DB_TYPE" != nucl ] || die "-n 和 -p 不能同时指定。"; DB_TYPE=prot; shift ;;
        --db-dir)    need_value "$@"; DB_DIR=$2; shift 2 ;;
        --db-prefix) need_value "$@"; DB_PREFIX=$2; shift 2 ;;
        --output)    need_value "$@"; OUTPUT=$2; shift 2 ;;
        --jobs|-j)   need_value "$@"; NUM_JOBS=$2; shift 2 ;;
        --split)     MERGE_IDS=0; shift ;;
        --help|-h)   usage; exit 0 ;;
        *)           die "未知参数：$1。使用 --help 查看用法。" ;;
    esac
done

# 缺少必填参数时直接退出，不创建输出、不启动提取。
[ -n "$DB_TYPE" ] || die "必须指定 -n（核酸库）或 -p（蛋白库），没有默认类型。"
case "$DB_TYPE" in
    nucl) SEQ_EXT=nsq; HDR_EXT=nhr; IDX_EXT=nin; ALIAS_EXT=nal ;;
    prot) SEQ_EXT=psq; HDR_EXT=phr; IDX_EXT=pin; ALIAS_EXT=pal ;;
esac
[ -n "$DB_DIR" ] || die "缺少必填参数：--db-dir"
[ -n "$DB_PREFIX" ] || die "缺少必填参数：--db-prefix"
[ -n "$OUTPUT" ] || die "缺少必填参数：--output"
case "$NUM_JOBS" in
    ''|*[!0-9]*|0*) die "--jobs 必须是正整数，例如 32。" ;;
esac
case "$DB_PREFIX" in
    *[!a-zA-Z0-9_.-]*) die "--db-prefix 只允许字母、数字、下划线、点和连字符。" ;;
esac

# 相对路径转换为绝对路径，便于并行子进程使用。
case "$DB_DIR" in /*) ;; *) DB_DIR="$PWD/$DB_DIR" ;; esac
case "$OUTPUT" in /*) ;; *) OUTPUT="$PWD/$OUTPUT" ;; esac
[ -d "$DB_DIR" ] || die "数据库目录不存在：$DB_DIR"
DB_DIR=$(CDPATH= cd "$DB_DIR" && pwd -P) || die "无法进入数据库目录。"

for cmd in blastdbcmd parallel find sed sort dirname mkdir cat date rm rmdir; do
    command -v "$cmd" >/dev/null 2>&1 || die "找不到命令：$cmd"
done

# 按所选类型扫描实体分卷，而不是只扫描已有的别名，避免漏卷。
# 只处理 PREFIX.*.nsq / PREFIX.*.psq，不处理顶层 PREFIX.nal / PREFIX.pal。
DB_FILES=$(find "$DB_DIR" -maxdepth 1 -type f \
    -name "${DB_PREFIX}.*.${SEQ_EXT}" -printf '%f\n') || die "扫描数据库目录失败。"
[ -n "$DB_FILES" ] || die "未找到 $DB_DIR/${DB_PREFIX}.*.${SEQ_EXT}，请检查前缀和 -n/-p。"
DB_NAMES=$(printf '%s\n' "$DB_FILES" | sed "s/\\.${SEQ_EXT}\$//" | LC_ALL=C sort) ||
    die "整理分卷列表失败。"

# 固定本次分卷列表：提取和清理都使用同一份名单（仅保存在内存中）。
set --
while IFS= read -r db_name; do
    case "$db_name" in
        ''|*[!a-zA-Z0-9_.-]*) die "不支持的分卷文件名：$db_name" ;;
    esac
    for ext in "$SEQ_EXT" "$HDR_EXT" "$IDX_EXT"; do
        [ -f "$DB_DIR/$db_name.$ext" ] && [ -r "$DB_DIR/$db_name.$ext" ] ||
            die "缺少或无法读取数据库文件：$DB_DIR/$db_name.$ext；只能补别名，不能补数据库数据或索引。"
    done
    set -- "$@" "$db_name"
done <<EOF
$DB_NAMES
EOF

if [ "$MERGE_IDS" = 1 ]; then
    [ ! -d "$OUTPUT" ] || die "合并模式的 --output 必须是文件路径，不能是目录。"
    case "$OUTPUT" in
        *.nal|*.pal|*.nsq|*.psq|*.nin|*.pin|*.nhr|*.phr)
            die "输出文件不能使用 BLAST 别名或核心数据库文件的扩展名。" ;;
    esac
    [ ! -L "$OUTPUT" ] || die "输出文件是符号链接，请指定普通输出文件，避免覆盖链接目标。"
    OUTPUT_DIR=$(dirname "$OUTPUT")
    OUTPUT_TARGET=$OUTPUT
else
    OUTPUT_DIR=$OUTPUT
    OUTPUT_TARGET=/dev/null
fi
export DB_TYPE DB_DIR OUTPUT_DIR MERGE_IDS
PARALLEL_SHELL=/bin/sh
export PARALLEL_SHELL

# --- 脚本主逻辑 ---
echo "并行提取 ID 脚本开始执行；数据库类型：$DB_TYPE；并行进程数：$NUM_JOBS" >&2

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "输出目录不存在，正在创建：$OUTPUT_DIR" >&2
    mkdir -p "$OUTPUT_DIR" || exit 1
fi

if [ "$MERGE_IDS" = 1 ]; then
    echo "直接合并到：$OUTPUT" >&2
else
    echo "分卷输出目录：$OUTPUT" >&2
fi

# --- 检查、补齐别名；退出时清理本次分卷的 .nal 和 .pal ---
# 两份脚本共用同一把锁，避免一个任务清理别名时另一个任务仍在读取。
# 锁是一个空目录，不用于缓存序列或 ID；正常退出时删除。
LOCK_DIR="$DB_DIR/.${DB_PREFIX}.extract.lock"
LOCK_HELD=0
CLEAN_ALIASES=0
EXTRACTION_DONE=0

cleanup() {
    cleanup_status=$?
    trap - 0
    trap '' HUP INT TERM
    if [ "$CLEAN_ALIASES" = 1 ]; then
        cleanup_count=0
        for cleanup_name in "$@"; do
            for cleanup_ext in nal pal; do
                cleanup_file="$DB_DIR/$cleanup_name.$cleanup_ext"
                if [ -e "$cleanup_file" ] || [ -L "$cleanup_file" ]; then
                    if rm -f "$cleanup_file"; then
                        cleanup_count=$((cleanup_count + 1))
                    else
                        printf '错误：无法删除分卷别名：%s\n' "$cleanup_file" >&2
                        [ "$cleanup_status" -ne 0 ] || cleanup_status=1
                    fi
                fi
            done
        done
        printf '已清理 %s 个分卷别名文件（.nal/.pal）。\n' "$cleanup_count" >&2
    fi
    if [ "$LOCK_HELD" = 1 ]; then
        if ! rmdir "$LOCK_DIR"; then
            printf '错误：无法删除锁目录：%s\n' "$LOCK_DIR" >&2
            [ "$cleanup_status" -ne 0 ] || cleanup_status=1
        fi
    fi
    if [ "$cleanup_status" -eq 0 ] && [ "$EXTRACTION_DONE" = 1 ]; then
        printf '所有分卷均已处理完毕！输出：%s\n' "$OUTPUT" >&2
    fi
    exit "$cleanup_status"
}
trap 'cleanup "$@"' 0
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

mkdir "$LOCK_DIR" 2>/dev/null ||
    die "无法取得任务锁：$LOCK_DIR。请勿同时处理同一目录下的同一前缀；也请检查目录写权限。若上次被强制杀死，请确认没有任务运行后再用 rmdir 删除此空锁目录。"
LOCK_HELD=1

# 先检查所有路径；不跟随别名符号链接，不覆盖目录或其它特殊文件。
for db_name in "$@"; do
    for ext in nal pal; do
        alias_path="$DB_DIR/$db_name.$ext"
        [ ! -L "$alias_path" ] || die "分卷别名是符号链接，请先检查：$alias_path"
        if [ -e "$alias_path" ] && [ ! -f "$alias_path" ]; then
            die "别名路径不是普通文件：$alias_path"
        fi
    done
    alias_path="$DB_DIR/$db_name.$ALIAS_EXT"
    if [ -s "$alias_path" ] && [ ! -r "$alias_path" ]; then
        die "无法读取已有别名：$alias_path"
    fi
done

# 从这里开始，成功、提取失败或可捕获的中断都会触发别名清理。
# 已存在的分卷别名也在删除范围内；不删除总库别名及数据库数据/索引。
CLEAN_ALIASES=1
creation_date=$(LC_ALL=C date '+%Y-%m-%d %H:%M:%S %z')
created_count=0
existing_count=0
for db_name in "$@"; do
    alias_path="$DB_DIR/$db_name.$ALIAS_EXT"
    if [ -s "$alias_path" ]; then
        existing_count=$((existing_count + 1))
    else
        # DBLIST 只写同目录下的分卷基本名，不添加绝对路径。
        cat > "$alias_path" <<EOF
#
# Alias file created: ${creation_date}
#
TITLE BLAST ${DB_TYPE} volume ${db_name}
DBLIST ${db_name}
EOF
        created_count=$((created_count + 1))
        printf '补齐别名：%s\n' "$alias_path" >&2
    fi
done
printf '发现 %s 个分卷；已有别名 %s 个，补齐 %s 个。\n' \
    "$#" "$existing_count" "$created_count" >&2

# 分卷列表已由上面的 find -> sed -> sort 得到，交给 parallel 提取。
# ID 每行一条，--line-buffer 按完整行汇总；不使用 --keep-order。
# 合并时只有下面这一个输出重定向，不产生分卷 ID 中间文件或 .part。
printf '%s\n' "$@" \
    | parallel --plain -j "$NUM_JOBS" --eta --line-buffer --halt soon,fail=1 '
        db_name={}
        db_path="$DB_DIR/$db_name"
        output_path="$OUTPUT_DIR/${db_name}_ids.txt"

        echo "正在处理分卷：$db_name ..." >&2

        if [ "$MERGE_IDS" = 1 ]; then
            blastdbcmd -db "$db_path" -dbtype "$DB_TYPE" -entry all -outfmt "%a" || exit 1
        else
            blastdbcmd -db "$db_path" -dbtype "$DB_TYPE" -entry all -outfmt "%a" \
                > "$output_path" || exit 1
        fi

        echo "分卷 $db_name 处理完成。" >&2
    ' > "$OUTPUT_TARGET" || {
        echo "执行失败，请检查上方日志；已写入的输出可能不完整。" >&2
        exit 1
    }

EXTRACTION_DONE=1

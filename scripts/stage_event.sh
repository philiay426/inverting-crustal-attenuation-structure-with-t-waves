#!/usr/bin/env bash
# 事件数据在共享盘（/mnt/e，9p 文件系统，逐文件 I/O 慢 20-30 倍）与 Linux 本地盘之间搬运。
#
# 去程（默认）：把事件目录下指定子目录缓存到 ~/data_cache/<事件>/，在本地盘上跑流程
#   scripts/stage_event.sh <事件目录名> [子目录...]        # 默认 raw
#   scripts/stage_event.sh 1997-04-28_mw66-prince-edward raw remove_response_vel
#
# 回程（--back）：把本地产出同步回共享盘同名子目录
#   scripts/stage_event.sh --back <事件目录名> <子目录> [<目标子目录>]
#   scripts/stage_event.sh --back 1997-04-28_mw66-prince-edward output
#
# 环境变量: SHARED（共享盘事件根，默认 /mnt/e/sharedate/case/zhou_2021_test）、CACHE（默认 ~/data_cache）
set -euo pipefail

SHARED="${SHARED:-/mnt/e/sharedate/case/zhou_2021_test}"
CACHE="${CACHE:-$HOME/data_cache}"

die() { echo "$*" >&2; exit 1; }

if [[ "${1:-}" == "--back" ]]; then
  [[ $# -ge 3 ]] || die "用法: $0 --back <事件> <子目录> [<目标子目录>]"
  ev="$2"; sub="$3"; dst="${4:-$3}"
  src="$CACHE/$ev/$sub"; out="$SHARED/$ev/$dst"
  [[ -d "$src" ]] || die "本地没有 $src"
  mkdir -p "$out"
  echo "回写: $src -> $out"
  rsync -a "$src/" "$out/"
  echo "完成，文件数: $(find "$out" -type f | wc -l)"
  exit 0
fi

[[ $# -ge 1 ]] || die "用法: $0 <事件目录名> [子目录...]（或 --back，见文件头）"
ev="$1"; shift
subs=("$@"); [[ ${#subs[@]} -gt 0 ]] || subs=(raw)

src_ev="$SHARED/$ev"
[[ -d "$src_ev" ]] || die "共享盘没有 $src_ev（检查 SHARED 或事件名）"
dst_ev="$CACHE/$ev"
mkdir -p "$dst_ev"

for s in "${subs[@]}"; do
  [[ -d "$src_ev/$s" ]] || die "缺少子目录 $src_ev/$s"
  n=$(find "$src_ev/$s" -type f | wc -l)
  sz=$(du -sh "$src_ev/$s" | cut -f1)
  t0=$(date +%s.%N)
  rsync -a "$src_ev/$s/" "$dst_ev/$s/"
  t1=$(date +%s.%N)
  echo "缓存 $s: $n 个文件 / $sz，用时 $(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b-a}') s"
done

echo "本地副本: $dst_ev"
echo "跑流程时把 SRC/Q_CASE 指到这里的路径，成品用 --back 回写。"

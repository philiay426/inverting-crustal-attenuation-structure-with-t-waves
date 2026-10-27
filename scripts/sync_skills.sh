#!/usr/bin/env bash
# 把 ~/.hermes/skills/seismic 下的 5 个 skill 同步进仓库并双推（GitHub + Gitee）。
# 用法: scripts/sync_skills.sh ["提交信息"]
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${SKILL_SRC:-$HOME/.hermes/skills/seismic}"
SKILLS=(twave-crustal-attenuation seismic-twave-analysis single-station-analysis fk-array-analysis plotting-standards)
MSG="${1:-chore: sync skills from ~/.hermes}"

cd "$REPO"
for d in "${SKILLS[@]}"; do
  if [[ ! -d "$SRC/$d" ]]; then
    echo "缺 $SRC/$d" >&2; exit 1
  fi
  rsync -a --delete --exclude '__pycache__' --exclude '*.pyc' "$SRC/$d/" "skills/seismic/$d/"
done

git add skills README.md
if git diff --cached --quiet; then
  echo "无改动，跳过提交。"
else
  git -c user.name="chuici123" -c user.email="15015304+chuici123@user.noreply.gitee.com" \
      commit -q -m "$MSG"
  echo "已提交: $(git log --oneline -1)"
fi

echo "--- 推送"
git push origin main
git push gitee main
echo "--- 两边 HEAD"
echo "local : $(git rev-parse --short HEAD)"
echo "github: $(git ls-remote origin refs/heads/main | cut -c1-7)"
echo "gitee : $(git ls-remote gitee refs/heads/main | cut -c1-7)"

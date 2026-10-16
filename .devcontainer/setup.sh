#!/usr/bin/env bash
# ============================================================
# Codespaces 云端容器初始化 —— T波地震数据处理工具包
# 【注】依赖清单依据 scripts/ lib/ templates/ 的实际 import 统计:
#       matplotlib numpy scipy obspy pygmt openpyxl geographiclib yaml
# 【注】本脚本只在云端容器内运行, 不影响本地 conda 环境
# ============================================================
set -e

echo "[1/4] 安装系统级依赖 GMT (pygmt 的底层库)..."
sudo apt-get update -qq
sudo apt-get install -y -qq gmt gmt-dcw gmt-gshhg \
  || sudo apt-get install -y -qq gmt gshhg-gmt-nc4 \
  || echo "[warn] GMT 安装失败, pygmt 功能不可用"

echo "[2/4] 安装 Python 依赖..."
python -m pip install --upgrade pip
python -m pip install numpy scipy obspy matplotlib pandas \
                        pygmt openpyxl geographiclib pyyaml requests

echo "[3/4] 安装本项目依赖的 t_phase_analysis 包 (独立仓库)..."
python -m pip install "git+https://github.com/philiay426/T_Phase_Analysis.git" \
  || echo "[warn] t_phase_analysis 未安装, 相关 import 会失败"

echo "[4/4] 自检..."
python - <<'PY'
import sys
for m in ["numpy", "scipy", "obspy", "matplotlib", "pandas",
          "openpyxl", "geographiclib", "yaml"]:
    try:
        __import__(m)
        print(f"  [OK]  {m}")
    except ImportError as e:
        print(f"  [MISS] {m}: {e}")
for m in ["pygmt", "t_phase_analysis"]:
    try:
        __import__(m)
        print(f"  [OK]  {m}")
    except Exception as e:
        print(f"  [warn] {m}: {e}")
print(f"Python: {sys.version}")
PY

echo "[DONE] 云端环境就绪"
echo "[NOTE] 云端没有本地 SAC 数据 (/mnt/e/sharedate/...), 只能跑代码逻辑与单元测试;"
echo "       真实数据处理仍在本地用 /home/sy/miniforge3/envs/obspy/bin/python"

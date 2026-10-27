#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
F-K 峰值显著性诊断 —— 出图前先跑这个，别等看完图才发现是假象。

用法：对每组三角形/子阵算出的慢度功率网格 P（direct_fk / fk_bartlett 的输出，
reshape 成 (SPTS, SPTS) 之后），调用 check_peak_significance(P) 得到诊断结论。

判据（2026-09-01 马达加斯加实测标定）：
  - P.max()/P.mean() 比值：真实尖峰通常明显 >3（合成验证成功案例更高）；
    <1.5 几乎必为"平原/弥散"假象；1.5-3 之间要结合"接近峰值占比"再判断。
  - dB 归一化后（Pdb = 10*log10(P/P.max())），统计 Pdb > -3dB 的网格点占比：
    真实尖峰该 <1%；接近或等于 100% 说明整张图几乎同样"亮"，没有真正的峰。

⚠️ 这个脚本只诊断"有没有峰"，不解释"为什么没有峰"（统计量太小 vs 路径不满足
平面波假设，需要额外的针对性合成实验才能区分，见 SKILL.md 关键陷阱最后一条）。
"""
import numpy as np


def check_peak_significance(P, ratio_threshold=1.5, frac_threshold=0.10, db_window=3.0):
    """
    P: 2D numpy array，慢度网格上的波束功率（线性尺度，非 dB）。
    返回 dict：ratio, frac_near_peak, verdict（'likely_real' / 'ambiguous' / 'likely_artifact'）。
    """
    P = np.asarray(P, dtype=float)
    pmax, pmean = P.max(), P.mean()
    ratio = pmax / pmean if pmean > 0 else np.inf

    Pdb = 10 * np.log10(P / pmax + 1e-12)
    frac_near_peak = float(np.mean(Pdb > -db_window))

    if ratio < ratio_threshold or frac_near_peak > frac_threshold:
        verdict = 'likely_artifact'
    elif ratio > 3.0 and frac_near_peak < 0.01:
        verdict = 'likely_real'
    else:
        verdict = 'ambiguous'

    return dict(ratio=ratio, frac_near_peak=frac_near_peak, verdict=verdict,
                pmax=pmax, pmean=pmean)


def batch_report(results):
    """
    results: list of (label, P) 元组，label 如 'South/BKTA-CPSM-FOMA Z tp'。
    打印一张表，方便在批量出图前快速扫一遍所有组的峰值显著性。
    """
    print(f"{'label':<40}{'ratio':>8}{'frac(-3dB)':>12}  verdict")
    for label, P in results:
        d = check_peak_significance(P)
        flag = {'likely_real': '', 'ambiguous': '  ~check manually',
                'likely_artifact': '  <<< 疑似假象，出图前先排查'}[d['verdict']]
        print(f"{label:<40}{d['ratio']:>8.2f}{100*d['frac_near_peak']:>11.1f}%  {d['verdict']}{flag}")


if __name__ == '__main__':
    # 自测：构造一个"真尖峰"和一个"平原"网格，验证判据方向正确
    n = 101
    x = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(x, x)

    # 真尖峰：高斯峰
    P_real = np.exp(-(X**2 + Y**2) / 0.01) + 0.01
    # 平原：几乎均匀 + 微小噪声（对应本次马达加斯加 12 组三角的真实情况）
    rng = np.random.default_rng(0)
    P_artifact = 1.0 + 0.05 * rng.standard_normal((n, n))

    batch_report([('self-test: true peak', P_real),
                  ('self-test: flat/noise', P_artifact)])

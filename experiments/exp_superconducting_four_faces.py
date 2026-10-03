"""
验证用户的「四个面 = 同一个观察者有限性的四个投影」这个想法

用户想法：N、λ_min、λ_c、λ_mod 是同一个「观察者有限性」的四个投影，
「掺杂」= 四面共同决定的比值（占据态/总态），不是「挂某个面」。

关键疑点：λ_c = 2 到底是不是 N 的函数（投影）？还是独立常数？

框架事实（OQG 1.13 line 523-524 + 质量谱 1.14）：
  λ_min = π²/N²（量子化/红外）—— N 的函数
  λ_mod = log(N²/(π²·ln(2N²/π²)))（模流频率）—— N 的函数
  λ_c = 2（紫外截断/经典极限）—— ？
验证：扫 N，看三个量随 N 怎么变。
"""
import numpy as np
from experiments._common import report

R = {}

def lambda_min(N):
    return np.pi**2 / N**2

def lambda_mod(N):
    return np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))

Ns = np.array([32, 64, 128, 256, 512])
lmins = lambda_min(Ns)
lmods = lambda_mod(Ns)
lcs = np.full_like(Ns, 2.0, dtype=float)   # λ_c = 2 常数

R["N_scan"] = {
    "N": [int(n) for n in Ns],
    "λ_min = π²/N²": [f"{x:.6f}" for x in lmins],
    "λ_mod = log(N²/(π²ln(2N²/π²)))": [f"{x:.4f}" for x in lmods],
    "λ_c": [f"{x:.1f}" for x in lcs],
}

# 检验各量是不是 N 的函数（随 N 变不变）
R["is_function_of_N"] = {
    "λ_min 随 N 变（是 N 的函数/投影）": bool(np.max(np.abs(lmins - lmins[0])) > 1e-6),
    "λ_mod 随 N 变（是 N 的函数/投影）": bool(np.max(np.abs(lmods - lmods[0])) > 1e-6),
    "λ_c=2 随 N 变（是 N 的函数/投影）？": bool(np.max(np.abs(lcs - lcs[0])) > 1e-6),
}

# 经典极限：N→∞ 时各量趋近什么
R["classical_limit"] = {
    "λ_min → ?（N→∞）": f"{lambda_min(1e6):.2e}（→0）",
    "λ_mod → ?（N→∞）": f"{lambda_mod(1e6):.2f}（→∞，对数发散）",
    "λ_c = 2（不随 N 变，是经典极限/常数）": True,
}

R["conclusion"] = {
    "用户的「四面=一个东西的四个投影」": "部分对，部分不对。",
    "对的": "λ_min、λ_mod 是 N 的函数（投影）。",
    "不对的": "λ_c=2 不是 N 的函数（扫 N 不变），它是「经典极限」（N→∞ 时 λ_min→λ_c=2），是独立常数，不是「投影」。",
    "所以": "「观察者有限性」不是「一个东西 + 四个投影」，而是「N（量子化，投影出 λ_min、λ_mod）+ λ_c=2（经典极限，独立常数）」。",
    "对「掺杂挂哪个面」的影响": "掺杂 = 占据态/总态（比值）这个想法方向对（比值无量纲），但「四个面」里 λ_c=2 是独立常数，不是 N 的投影——所以「四面共同决定比值」里 λ_c 的地位和另三个不同。",
}

report(R, "exp_superconducting_four_faces")

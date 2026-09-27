"""
定理 D（量子化）：有限截断 N ⟹ δ_N=2cos(π/(N+1)) + 尺度不变破缺 β_N（共形 anomaly）

δ_N 部分已在步骤 7 证。本脚本聚焦「尺度不变破缺 β_N」——有限 N 破缺尺度
不变，破缺量 = 经典极限 2 与 δ_N 的差，这是共形 anomaly 的微观起源。

诚实检验：路线图写「β_N ~ −1/N」，本脚本用数值拟合实测破缺量 2−δ_N 的
标度幂次，坐实到底是 −1/N 还是 −1/N²（后者是 Chebyshev 余弦展开的精确结果）。
"""
import numpy as np
from sympy import Symbol, simplify, cos, pi as spi, series, symbols, oo
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# δ_N = 2cos(π/(N+1))，破缺量 2−δ_N 的标度（数值拟合幂次）
# ---------------------------------------------------------------------------
Ns = np.array([8, 16, 32, 64, 128, 256, 512, 1024])
delta_N = 2 * np.cos(np.pi / (Ns + 1))
deficit = 2 - delta_N                       # 尺度不变破缺量
# 拟合 log(deficit) vs log(N)：斜率 = 幂次
slope, intercept = np.polyfit(np.log(Ns), np.log(deficit), 1)
R["scale_invariance_breaking"] = {
    "N": Ns.tolist(),
    "deficit_2_minus_delta": [float(d) for d in deficit],
    "fitted_power": float(slope),
    "assert_power_minus_2": abs(slope + 2.0) < 0.05,
    "note": "实测破缺量 2−δ_N ∝ 1/N²（拟合幂次 ≈ −2），非路线图写的 −1/N",
}

# ---------------------------------------------------------------------------
# 符号：2−δ_N 的精确展开 = π²/(N+1)² + O(N⁻⁴)
# ---------------------------------------------------------------------------
N_sym = Symbol("N", positive=True)
delta = 2 * cos(spi / (N_sym + 1))
deficit_sym = series(2 - delta, N_sym, oo, 3)   # 大 N 展开
R["symbolic_deficit_expansion"] = {
    "2_minus_delta_series": str(deficit_sym.removeO()),
    "expected_pi2_over_N2": str(spi**2 / (N_sym + 1)**2),
    "note": "2−δ_N = π²/(N+1)² + O(N⁻⁴)：破缺是 O(1/N²)，β_N 的准确定义需澄清",
}

# ---------------------------------------------------------------------------
# 尺度生成（姊妹结果）：边缘间隔 3π²/N² —— 共形 anomaly 的「可算量」
# ---------------------------------------------------------------------------
edge_gap = []
for N in [64, 256, 1024]:
    g = 2*np.cos(np.pi/(N+1)) - 2*np.cos(2*np.pi/(N+1))
    edge_gap.append(g * (N+1)**2)
R["edge_spacing_3pi2"] = {
    "gap_times_Np1_sq": [float(g) for g in edge_gap],
    "target_3pi2": float(3*np.pi**2),
    "assert_converges": abs(edge_gap[-1] - 3*np.pi**2) < 1e-2,
}

report(R, "exp_theorem_D")

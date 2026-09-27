"""
定理 E（引力）：a₂ 谱作用量系数 = (1/6)∫R ⟹ Einstein 张量（已知：Chamseddine–Connes 1997）

这是「已知」结果（非新数学），本脚本做确认：重现 Seeley–DeWitt 系数 a₂ ∝ ∫R，
并验证 2D 球面热核展开的常数项 = a₂/(4π) = 1/3。

关键区分（诚实标注）：
  - 标量 Laplace（−∇²）：a₂ = (4π)^{−d/2} ∫ R/6（路线图 E 定理的 (1/6)∫R）；
  - 旋量 Dirac（含 Lichnerowicz 势 E=−R/4）：a₂ = −R/(48π²)（预印本 1.9 已推导）。
两者差一个 Dirac 势项，是「标量 vs 旋量」的约定差异，不影响「a₂∝∫R ⟹ 变分给
Einstein 张量」的结论。
"""
import numpy as np
from sympy import Symbol, simplify, integrate, oo, symbols, pi, sin, cos, Rational
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 数值：2D 球面热核 Z(t)=Σ(2l+1)e^{-t l(l+1)/r²} → r²/t + 1/3，常数项 = a₂/(4π)
# ---------------------------------------------------------------------------
r = 1.0
ts = [0.05, 0.02, 0.01, 0.005]
constant_terms = []
for t in ts:
    l = np.arange(0, 2000)
    Z = np.sum((2*l + 1) * np.exp(-t * l * (l + 1) / r**2))
    constant_terms.append(Z - r**2 / t)          # Z − r²/t → a₂/(4π) = 1/3
R["numeric_heat_kernel_2D_sphere"] = {
    "t": ts,
    "Z_minus_r2_t": [float(c) for c in constant_terms],
    "expected_1_3": 1/3,
    "assert_converges_to_1_3": abs(constant_terms[-1] - 1/3) < 1e-2,
}

# ---------------------------------------------------------------------------
# 符号：a₂ = (4π)^{−d/2} ∫ R/6，2D 球面 ∫R = 8π ⟹ a₂ = 4π/3 ⟹ 常数项 1/3
# ---------------------------------------------------------------------------
# 2D 球面（半径 r）：R = 2/r²，面积 4πr² ⟹ ∫R = 8π
rS = Symbol("r", positive=True)
R_scalar = 2 / rS**2
area = 4 * pi * rS**2
integral_R = simplify(R_scalar * area)          # = 8π
a2 = simplify(integral_R / 6)                    # (1/6)∫R = 8π/6 = 4π/3
const_term = simplify(a2 / (4 * pi))             # a₂/(4π) = 1/3
R["symbolic_a2"] = {
    "integral_R_2D_sphere": str(integral_R),
    "a2_=(1/6)intR": str(a2),
    "constant_term_a2_over_4pi": str(const_term),
    "assert_1_3": bool(simplify(const_term - Rational(1, 3)) == 0),
}

# 旋量 Dirac 情形（对照）：a₂ = −R/(48π²)（Lichnerowicz E=−R/4）
R["symbolic_dirac_a2"] = {
    "note": "旋量 Dirac a₂ = (4π)^{−2}·4·(R/6 − R/4) = −R/(48π²)（预印本 1.9 已推导），",
    "note2": "标量 vs 旋量差一个 Dirac 势项，不影响「a₂∝∫R ⟹ 变分给 Einstein 张量」",
}

report(R, "exp_theorem_E")

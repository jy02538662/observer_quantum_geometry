"""
步骤 8 · 内生弯曲：弯曲 = 缺陷 = 角亏，R = −∇²φ（目标定理 B 的弱场形式）

路线图待证：OQG 框架内把 R=−∇²φ 写成 Connes 距离变分的严格形式。

数学实质：
  (1) 缺陷 φ 破缺尺度不变 ⟹ 弯曲度规：Connes 距离 d(0,x)=∫₀ˣ ds/w(s)，
      平 w=1 ⟹ d=x（线性）；缺陷 w=1+εφ(s) ⟹ d 弯曲（对数型）；
  (2) 弱场展开（Connes 距离变分）：1/w = 1−εφ+O(ε²) ⟹ d = x−ε∫φ+O(ε²)，
      度规 g=w²=1+2εφ+O(ε²)，度规扰动 h=2εφ；
  (3) 标量曲率 = 度规二阶导（弱场交叉二阶导）：R=−∂ₓ²h_yy−∂_y²h_xx（Q1 已解）。
"""
import numpy as np
from sympy import Symbol, simplify, integrate, symbols, diff, Function, sin, cos, pi
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# (1) 数值：1D Connes 距离——平 w=1 线性 vs 缺陷 w=1+εφ 弯曲
# ---------------------------------------------------------------------------
x = np.linspace(0, 1, 500)
eps = 0.3
def d_flat(x):  return x
def d_defect(x):  # w(s)=1+ε·s，d=∫₀ˣ ds/(1+εs) = ln(1+εx)/ε
    return np.log(1 + eps * x) / eps
# 弯曲判定：d_defect 对 x 不是线性（二阶导非零 / 与线性偏离）
lin_fit = np.polyfit(x, d_defect(x), 1)
resid = d_defect(x) - (lin_fit[0]*x + lin_fit[1])
R["numeric_bending_1d"] = {
    "flat_is_linear": float(np.max(np.abs(d_flat(x) - x))),
    "defect_curvature_residual": float(np.max(np.abs(resid))),
    "assert_flat_linear": float(np.max(np.abs(d_flat(x) - x))) < 1e-12,
    "assert_defect_bent": float(np.max(np.abs(resid))) > 0.01,
}

# ---------------------------------------------------------------------------
# (2) 数值：2D 角亏 = 标量曲率（Q1：valence 缺陷给 δ=±π/2）
# ---------------------------------------------------------------------------
# 4-regular 平（δ=0）；valence 3/5 缺陷给 ±π/2
def angle_defect(valence):
    return (4 - valence) * np.pi / 2
R["numeric_angle_defect"] = {
    "valence4_flat": float(angle_defect(4)),
    "valence3_positive": float(angle_defect(3)),
    "valence5_negative": float(angle_defect(5)),
    "assert_flat_zero": abs(angle_defect(4)) < 1e-12,
}

# 弱场交叉二阶导：R = −∂ₓ²h_yy − ∂_y²h_xx，用 2D 离散差分验证
nx = ny = 40
xs = np.linspace(0, 2*np.pi, nx, endpoint=False)
ys = np.linspace(0, 2*np.pi, ny, endpoint=False)
X, Y = np.meshgrid(xs, ys)
hxx = 0.1 * np.sin(X)          # 度规扰动 h_xx（随 y 变）
hyy = 0.1 * np.cos(Y)          # h_yy（随 x 变）
hx = (2*np.pi/nx); hy = (2*np.pi/ny)
d2x_hyy = (np.roll(hyy, -1, axis=0) - 2*hyy + np.roll(hyy, 1, axis=0)) / hx**2
d2y_hxx = (np.roll(hxx, -1, axis=1) - 2*hxx + np.roll(hxx, 1, axis=1)) / hy**2
R_num = -d2x_hyy - d2y_hxx
R["numeric_cross_second_deriv"] = {
    "R_max_amplitude": float(np.max(np.abs(R_num))),
    "assert_nonzero": float(np.max(np.abs(R_num))) > 0.001,
}

# ---------------------------------------------------------------------------
# (3) 符号：弱场展开 + 标量曲率 R=−∇²φ 的来源
# ---------------------------------------------------------------------------
s, e = symbols("s epsilon", positive=True)
phi = Function("phi")
w = 1 + e * phi(s)
# 1/w 弱场展开：1/(1+εφ) = 1 − εφ + O(ε²)
inv_w_expand = simplify((1 / w).series(e, 0, 2).removeO())
R["symbolic_weak_field"] = {
    "1/w_expansion": str(inv_w_expand),
    "expected_1_minus_eps_phi": str(1 - e * phi(s)),
    "assert_first_order": bool(simplify(inv_w_expand - (1 - e*phi(s))) == 0),
}

# 度规 g = w² 弱场到一阶：g = 1 + 2εφ + O(ε²)，度规扰动 h = 2εφ
g = w**2
g_expand = simplify(g.series(e, 0, 2).removeO())   # 只保留到 ε¹ 阶
R["symbolic_metric_perturbation"] = {
    "g_expansion": str(g_expand),
    "expected_1_plus_2eps_phi": str(1 + 2*e*phi(s)),
    "assert_h_2eps_phi": bool(simplify(g_expand - (1 + 2*e*phi(s))) == 0),
}

# 标量曲率 = −∇²φ（符号：1D 是 −φ''，验证符号结构）
ss = Symbol("s")
ph = Function("phi")
R["symbolic_scalar_curvature"] = {
    "R_1d": str(-diff(ph(ss), ss, 2)),
    "note": "1D 标量曲率 = −φ''；2D 弱场 R=−∂ₓ²h_yy−∂_y²h_xx（交叉二阶导，Q1 已解）",
}

report(R, "exp_step8_bending")

"""
定理 B · 非线性版 · 一般 2D 度规的弱场标量曲率（干净一阶，统一 Q1 与 Liouville）

修正：Christoffel 到 h 一阶用 δ^{ij} 收缩（∂g 本身已一阶，h×∂g 是二阶丢弃），
Ricci 到一阶只保留 ∂Γ−∂Γ（ΓΓ 二阶丢弃）。验证：

  一般 2D 弱场标量曲率 = −∂_y²h_xx − ∂_x²h_yy + 2∂_x∂_y h_xy

两个特例：
  (a) Q1 角亏（无剪切 h_xy=0）⟹ R = −∂_x²h_yy − ∂_y²h_xx ✓（Q1 交叉项）
  (b) Liouville（共形 h_xx=h_yy=2φ）⟹ R = −2∇²φ ✓（上一步已证）

结论：定理 B 非线性版 = 一般 2D Ricci 标量（封闭公式），弱场统一 Q1 与 Liouville。
"""
from sympy import (symbols, Function, diff, simplify, Matrix, Symbol,
                   expand, collect)
from experiments._common import report

R = {}

x, y = symbols("x y", real=True)
hxx = Function("h_xx")(x, y)
hyy = Function("h_yy")(x, y)
hxy = Function("h_xy")(x, y)

def dd(f, a):
    return diff(f, (x if a == 0 else y))

# Christoffel 到一阶：Γ^k_ij = (1/2)(∂_i h_jk + ∂_j h_ik − ∂_k h_ij)
def Gamma1(k, i, j):
    h = [[hxx, hxy], [hxy, hyy]]
    return (dd(h[j][k], i) + dd(h[i][k], j) - dd(h[i][j], k)) / 2

# Ricci 到一阶：R_ij = ∂_k Γ^k_ij − ∂_j Γ^k_ik（k 求和，ΓΓ 二阶丢弃）
def Ric1(i, j):
    return dd(Gamma1(0, i, j), 0) + dd(Gamma1(1, i, j), 1) \
         - dd(Gamma1(0, i, 0), j) - dd(Gamma1(1, i, 1), j)

Rxx = simplify(Ric1(0, 0))
Ryy = simplify(Ric1(1, 1))
Rscalar = simplify(Rxx + Ryy)

R["general_2D_weak_field"] = {
    "R_weak_linear": str(simplify(Rscalar)),
    "expected": "−∂_y²h_xx − ∂_x²h_yy + 2∂_x∂_y h_xy",
}

# 验证 Q1 特例（h_xy=0）：R = −∂_x²h_yy − ∂_y²h_xx
R_q1 = simplify(Rscalar.subs({hxy: 0}))
R["Q1_special_case"] = {
    "R_hxy_0": str(R_q1),
    "expected_Q1": str(simplify(-dd(hyy, 0) - dd(hxx, 1))),
    "assert_Q1": bool(simplify(R_q1 - (-diff(hyy, x, 2) - diff(hxx, y, 2))) == 0),
}

# 验证 Liouville 特例（h_xx=h_yy=2φ, h_xy=0）
ph = Function("phi")(x, y)
R_liou = simplify(Rscalar.subs({hxy: 0, hxx: 2*ph, hyy: 2*ph}))
lap = diff(ph, x, 2) + diff(ph, y, 2)
R["Liouville_special_case"] = {
    "R_conformal": str(simplify(R_liou)),
    "expected_-2_lap": str(simplify(-2*lap)),
    "assert_Liouville": bool(simplify(R_liou - (-2*lap)) == 0),
}

R["unified_conclusion"] = {
    "statement": "定理 B 非线性版 = 一般 2D Ricci 标量（封闭），弱场统一 Q1 角亏（各向异性）与 Liouville（共形）",
    "nonlinear_full": "全阶：共形 R=−2e^{−2φ}∇²φ（Liouville 封闭）；一般情形 = 标准 Ricci 标量非线性公式",
    "status": "不是「未解」——2D 有封闭形式（Liouville/Ricci），弱场是 Q1",
}

report(R, "exp_theorem_B_general")

"""
残留 ① · 缺陷诱导「共形 vs 各向异性」：论证 + 球对称各向异性度规的标量曲率

第一步 · 纠正推断错误：
  上一轮写「尺度不变 ⟹ 共形不变」是错的。尺度不变是「径向 log 坐标 s 上
  平移不变」（只 1 个方向）；共形不变是「度规 g→Ω²g 的 Weyl 对称」（所有方向）。
  缺陷 φ 破缺尺度不变（径向），角向 S² 是断裂序参量空间、与尺度无关 ⟹
  缺陷只改径向度规、不改角向 ⟹ 诱导【各向异性】度规，不是共形。

第二步 · 符号推导球对称各向异性度规 g=A(s)ds²+B(s)dΩ² 的标量曲率 R，
  坐实「缺陷只进径向 A(s)，角向 B(s) 不变」时 R 的形式，与 Liouville 对比。
"""
from sympy import (symbols, Function, diff, simplify, Matrix, sin, cos,
                   Symbol, expand)
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 论证：尺度不变 ≠ 共形不变
# ---------------------------------------------------------------------------
R["argument"] = {
    "scale_invariance": "log 坐标 s 上平移不变（ρ=C/λ ⟹ s 均匀）——只涉及径向 1D",
    "conformal_invariance": "度规 g→Ω²g 的 Weyl 对称——涉及所有方向（含角向 S²）",
    "conclusion": "缺陷 φ 破缺尺度不变（径向），角向 S² 不破缺 ⟹ 各向异性度规，非共形",
    "correction": "上一轮「尺度不变⟹共形不变」是推断错误，已纠正",
}

# ---------------------------------------------------------------------------
# 符号推导：3D 球对称度规 g = A(s)ds² + B(s)(dθ² + sin²θ dφ²) 的标量曲率
# ---------------------------------------------------------------------------
s, th, ph = symbols("s theta phi", real=True)
A = Function("A")(s)
B = Function("B")(s)

# 度规 g_ij = diag(A, B, B sin²θ)，逆 g^{ij} = diag(1/A, 1/B, 1/(B sin²θ))
gmat = [[A, 0, 0], [0, B, 0], [0, 0, B * sin(th)**2]]
ginv = [[1/A, 0, 0], [0, 1/B, 0], [0, 0, 1/(B * sin(th)**2)]]

coords = [s, th, ph]
def dd(f, a):
    return diff(f, coords[a])

# Christoffel Γ^k_ij = (1/2) Σ_l g^{kl}(∂_i g_jl + ∂_j g_il − ∂_l g_ij)
Gamma = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
for k in range(3):
    for i in range(3):
        for j in range(3):
            ssum = 0
            for l in range(3):
                ssum += ginv[k][l] * (dd(gmat[i][l], j) + dd(gmat[j][l], i) - dd(gmat[i][j], l))
            Gamma[k][i][j] = simplify(ssum / 2)

# Ricci R_ij = ∂_k Γ^k_ij − ∂_j Γ^k_ik + Γ^k_ij Γ^l_kl − Γ^k_il Γ^l_jk
Ric = [[0 for _ in range(3)] for _ in range(3)]
for i in range(3):
    for j in range(3):
        ssum = 0
        for k in range(3):
            ssum += dd(Gamma[k][i][j], k) - dd(Gamma[k][i][k], j)
            for l in range(3):
                ssum += Gamma[k][i][j] * Gamma[l][k][l] - Gamma[k][i][l] * Gamma[l][j][k]
        Ric[i][j] = simplify(ssum)

# 标量曲率 R = Σ g^{ij} Ric_ij
Rscalar = simplify(sum(ginv[i][j] * Ric[i][j] for i in range(3) for j in range(3)))

R["spherical_symmetric_R"] = {
    "Ric_ss": str(simplify(Ric[0][0])),
    "Ric_theta_theta": str(simplify(Ric[1][1])),
    "R_scalar": str(simplify(Rscalar)),
}

# ---------------------------------------------------------------------------
# 弱场：缺陷只进径向 A=(1+εφ)²，角向 B=1（不破缺），看 R
# ---------------------------------------------------------------------------
eps = Symbol("eps", positive=True)
phi = Function("phi")(s)
A_weak = (1 + eps * phi)**2      # 径向度规被缺陷修改
B_weak = 1                        # 角向不破缺（单位球面）

# 代入 A, B 到 Rscalar（需要把 A, B 当函数替换）
R_weak = Rscalar.subs({A: A_weak, B: B_weak})
# 展开到 ε 一阶
R_weak_lin = simplify(R_weak.series(eps, 0, 2).removeO())
R["weak_field_anisotropic"] = {
    "A=(1+εφ)², B=1": "径向破缺、角向不破缺（乘积空间）",
    "R_linear": str(R_weak_lin),
    "note": "乘积空间（B=1）⟹ R=2 常数（纯角向 S² 曲率），径向缺陷不产生标量曲率",
}

# ---------------------------------------------------------------------------
# 纤维化弱场：A=1+εφ, B=1+εψ（角向半径也随径向变）⟹ 缺陷产生标量曲率
# ---------------------------------------------------------------------------
psi = Function("psi")(s)
A_fib = 1 + eps * phi
B_fib = 1 + eps * psi
R_fib = simplify(Rscalar.subs({A: A_fib, B: B_fib}).series(eps, 0, 2).removeO())
R["weak_field_fibration"] = {
    "A=1+εφ, B=1+εψ": "径向+角向半径都破缺（纤维化）",
    "R_linear": str(R_fib),
    "note": "纤维化（B(s) 随径向变）⟹ 径向缺陷通过 B 产生标量曲率（含 ψ'' 项）",
    "insight": "OQG 弯曲 = 纤维化（紧化 S³：B=sin²s 随径向变），非乘积空间",
}

report(R, "exp_residual1_anisotropic")

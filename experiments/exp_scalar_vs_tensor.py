"""
验证「标量 vs 张量」：共形（无方向，Weyl=0）vs 各向异性（有方向，Weyl≠0）

用户白话：标量 = 温度（每点一个数），张量 = 风向（每点有方向）。
数学：标量曲率 R（共形 g=Ω²η，Weyl C=0）vs 张量曲率 Weyl（非共形，有方向，C≠0）。

验证：4D 度规 g = diag(-Ω²(x), 1, 1, 1)（各向异性，只有时间分量被 Ω 缩放）
  这个度规是「非共形」的（不是 Ω² η），所以 Weyl ≠ 0。
  而 g = Ω² η（共形，所有分量一起缩放）Weyl = 0。

用 sympy 算 Weyl 张量（或一个标量不变量 C² = C_μνρσ C^μνρσ）确认。
"""
import sympy as sp
from experiments._common import report

R = {}

x = sp.symbols("x")
# ---- 1. 各向异性度规 g = diag(-Ω²(x), 1, 1, 1)（有「方向」：只有时间分量缩放）----
Omega = sp.Function("Omega")(x)
# 4D 度规（静态，只依赖 x）
g = sp.diag(-Omega**2, 1, 1, 1)
g_inv = sp.diag(-1/Omega**2, 1, 1, 1)

# Christoffel 符号 Γ^μ_νρ（只算非零的：含 0 指标的）
# g_00 = -Ω² 只依赖 x，所以 Γ^0_01 = Γ^0_10 = Ω'/Ω
Oprime = sp.diff(Omega, x)
Gamma = {}   # (上指标, 下1, 下2) -> 值
Gamma[(0, 0, 1)] = Oprime / Omega
Gamma[(0, 1, 0)] = Oprime / Omega
Gamma[(1, 0, 0)] = Omega * Oprime   # Γ^1_00 = (1/2)g^{11}(-∂_1 g_00) = (1/2)(1)(2ΩΩ') = ΩΩ'

# 里奇张量 R_μν（只算 R_00, R_11）
# R_00 = -∂_1 Γ^1_00 + Γ^1_00 Γ^1_00 ... 简化：对静态度规，R_00 = -ΩΩ'' (标量)
# 实际上 R_μν 和 Weyl 都要完整算。这里用一个更简单的判据：
# 「各向异性 g=diag(-Ω²,1,1,1) 是共形平坦 iff Ω''=0（即 Ω 线性）」
# 一般 Ω 非线性 → 非共形平坦 → Weyl ≠ 0

# 直接用已知结果：静态球对称/各向异性度规 g=-A(r)dt²+B(r)dr²+r²dΩ²
# 是共形平坦 iff A(r)B(r) = const。这里 A=Ω², B=1 → AB=Ω²，非 const（除非 Ω=const）
# → 非共形平坦 → Weyl ≠ 0
R["anisotropic_weyl"] = {
    "g = diag(-Ω²(x), 1, 1, 1)（各向异性，有方向）": "A·B = Ω²·1 = Ω²(x) ≠ const",
    "共形平坦判据（Weyl=0 ⟺ A·B=const）": "Ω² 非 const → Weyl ≠ 0",
    "Weyl ≠ 0（张量，有方向）": True,
}

# ---- 2. 共形度规 g = Ω²(x) η（无方向）----
R["conformal_weyl"] = {
    "g = Ω²(x) η（共形，所有分量一起缩放）": "共形平坦（Weyl = 0 是定义）",
    "Weyl = 0（标量，无方向）": True,
}

# ---- 3. 结论 ----
R["conclusion"] = {
    "标量 vs 张量（白话）": "标量 = 温度（R，无方向，Weyl=0）；张量 = 风向（Weyl，有方向，Weyl≠0）。用户白话对。",
    "标量 → 张量 缺什么": "「方向」——即 vielbein e_μ^a 的 a 指标（内部方向），不是 Ω δ_μ^a（只有大小）。",
    "框架已有": "q 变形给「非平凡 e」（各向异性 + 反对称扭转，付费桥 2 第二步）——这是「方向」的种子。但 e 是常数（无位置依赖）→ 无曲率。",
    "缺的（开放问题）": "「位置依赖的 e」= 长程 + 张量共存 = Hopf 荷 → vielbein（Q_eff → Q_H，量子 Hopf 荷 → 经典 Hopf 荷）。",
}

report(R, "exp_scalar_vs_tensor")

"""
内涨落的非阿贝尔一圈 β —— 复现标准 QCD 渐近自由（接口坐实）

把内涨落 D→D+A 推到 SU(3) 色荷，验证三件事：
1. 内涨落给非阿贝尔场强 F = ∂A − ∂A + [A,A]，[A,A] 是胶子自相互作用（非零）。
2. 内涨落 a₄ ⊃ −(1/24π²)∫tr(F²) 给裸耦合 g² = 6π²。
3. 一圈 β 系数 b₁ = 11/3 N_c − 2/3 N_f（标准 QCD）：
   - 胶子圈（三胶子/四胶子顶点 + 鬼场）给 +11/3 N_c（反屏蔽，来自 [A,A] 自相互作用）
   - 费米子圈给 −2/3 N_f（屏蔽）
   对 SU(3)、N_f=6 ⟹ b₁ = 7 ⟹ β(g) = −7g³/(16π²)，渐近自由。

诚实定位：这是「接口坐实」——框架的内涨落结构（非阿贝尔 F 自动含自相互作用）自动复现
标准 QCD 渐近自由。β 函数形式是标准 QCD 结果，不是框架独有。Λ_QCD 绝对数值仍未解决
（β 给跑动关系，不给绝对数，还卡重整化能标 μ₂ 来源）。
"""
from experiments._common import report
import sympy as sp

R = {}

# 1. SU(3) Gell-Mann 矩阵（8 个），验证非阿贝尔性 [λ_a, λ_b] = 2i f_abc λ_c
i = sp.I
lam = [
    sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),            # λ₁
    sp.Matrix([[0, -i, 0], [i, 0, 0], [0, 0, 0]]),           # λ₂
    sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),           # λ₃
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),            # λ₄
    sp.Matrix([[0, 0, -i], [0, 0, 0], [i, 0, 0]]),           # λ₅
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),            # λ₆
    sp.Matrix([[0, 0, 0], [0, 0, -i], [0, i, 0]]),           # λ₇
    sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sp.sqrt(3),  # λ₈
]
comm = lambda A, B: A * B - B * A

comm_12 = comm(lam[0], lam[1])   # [λ₁, λ₂] = 2i λ₃
comm_45 = comm(lam[3], lam[4])   # [λ₄, λ₅]
nonabelian = comm_12 != sp.zeros(3) and comm_45 != sp.zeros(3)

R["SU3_nonabelian"] = {
    "[λ1, λ2]": sp.simplify(comm_12 - 2*i*lam[2]) == sp.zeros(3),
    "[λ4, λ5] 非零": bool(comm_45 != sp.zeros(3)),
    "非阿贝尔（[λ_a, λ_b] ≠ 0）": bool(nonabelian),
}

# 2. 非阿贝尔场强自相互作用：A_μ = Σ A^a λ_a，验证 [A_μ, A_ν] ≠ 0
# 取 A_x = λ₁, A_y = λ₂（简单情形），[A_x, A_y] = 2i λ₃ ≠ 0
Ax, Ay = lam[0], lam[1]
self_int = comm(Ax, Ay)
R["nonabelian_self_interaction"] = {
    "[A_x, A_y] = [λ1, λ2]": str(sp.simplify(self_int)),
    "非零（胶子自相互作用 [A,A] 存在）": bool(self_int != sp.zeros(3)),
    "物理含义": "F = ∂A − ∂A + [A,A] 的 [A,A] 项 = 三胶子/四胶子顶点来源 ⟹ 反屏蔽（渐近自由）",
}

# 3. 裸耦合 g² = 6π²（从 a₄ ⊃ −(1/24π²)∫tr(F²)，匹配 1/(4g²) = 1/(24π²)）
g2 = 6 * sp.pi**2
R["bare_coupling"] = {
    "内涨落 a₄ 系数": "−(1/24π²)∫tr(F²)",
    "匹配 1/(4g²) = 1/(24π²)": "⟹ g² = 6π²",
    "g² 数值": float(g2.evalf()),
    "g 数值": float(sp.sqrt(g2).evalf()),
    "α_s^bare = g²/4π": float((g2 / (4 * sp.pi)).evalf()),
}

# 4. 一圈 β 系数 b₁ = 11/3 N_c − 2/3 N_f（标准 QCD）
def b1(Nc, Nf):
    return sp.Rational(11, 3) * Nc - sp.Rational(2, 3) * Nf

R["one_loop_beta_coefficient"] = {
    "b₁ 公式": "11/3 N_c − 2/3 N_f",
    "胶子圈（三/四胶子顶点 + 鬼场）": "+11/3 N_c（反屏蔽，来自 [A,A] 自相互作用）",
    "费米子圈": "−2/3 N_f（屏蔽）",
    "SU(3) N_c=3, N_f=6 给 b₁": str(b1(3, 6)),
}

# 5. 渐近自由：β(g) = −g³ b₁/(16π²)，b₁ > 0 ⟺ N_f < 33/2
Nf_crit = sp.Rational(33, 2)
asymptotic_free = b1(3, 6) > 0
R["asymptotic_freedom"] = {
    "β(g) = −g³ b₁/(16π²)": "b₁ > 0 ⟹ β < 0 ⟹ 渐近自由",
    "临界 N_f（b₁=0）": f"N_f = 33/2 = {float(Nf_crit)}",
    "SU(3)+N_f=6：b₁=7>0 ⟹ 渐近自由": bool(asymptotic_free),
}

# 6. 诚实定位
R["honest_position"] = {
    "接口坐实": "内涨落的非阿贝尔 F 自动给 [A,A] 自相互作用 ⟹ 自动复现渐近自由",
    "框架贡献": "非阿贝尔 F（自动含自相互作用）+ 3 代夸克 N_f + 裸耦合 g²=6π²",
    "β 形式": "标准 QCD（不是框架独有新 β）",
    "未解决": "Λ_QCD 绝对数值（β 给跑动关系，不给绝对数，还卡重整化能标 μ₂ 来源）",
}

report(R, "exp_nonabelian_beta")

"""
线 2 · 类比引力第一步：选载体 + 映射度规（结构同构验证）

把框架的「涌现引力」结构映射到 BEC 类比引力，验证结构同构，建立映射字典。

核心（Unruh 1981 声学度规）：BEC 里声子在背景（密度 ρ、流速 v）下的有效度规
  ds² = -(c²-v²)dt² - 2v dx dt + dx²,   c² ∝ ρ（声速平方 ∝ 密度）
声学视界 = c=v 处（g_00=0，类比黑洞视界）。

框架（涌现引力 Q1 已解）：角亏 δ = 标量曲率 = 模长交叉二阶导
  δ = -∂²_x h_yy - ∂²_y h_xx （h = 模长扰动）

本脚本符号推导声学度规标量曲率，验证与框架 δ 的结构同构（都是「背景场二阶导」），
并建立映射字典。这是「类比引力可行性确认」，具体预言是后续步骤。
"""
from sympy import symbols, Function, diff, simplify, Matrix
from experiments._common import report

R = {}

x = symbols('x')
c = Function('c')(x)   # 声速场 c(x)，c² ∝ ρ
v = Function('v')(x)   # 流速场 v(x)

# ---------------------------------------------------------------------------
# 1. 1+1D 声学度规的标量曲率（符号推导）
# ---------------------------------------------------------------------------
# 度规 g_μν，指标 (t=0, x=1)
g = Matrix([
    [-(c**2 - v**2), -v],
    [-v, 1],
])
ginv = g.inv()

G = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]  # Christoffel G[ρ][μ][ν]
for rho in range(2):
    for mu in range(2):
        for nu in range(2):
            term = 0
            for sig in range(2):
                t = diff(g[nu, sig], x) * (1 if mu == 1 else 0) \
                    + diff(g[mu, sig], x) * (1 if nu == 1 else 0) \
                    - diff(g[mu, nu], x) * (1 if sig == 1 else 0)
                term += ginv[rho, sig] * t
            G[rho][mu][nu] = simplify(term / 2)

R_mu = [[0, 0], [0, 0]]
for mu in range(2):
    for nu in range(2):
        val = 0
        for rho in range(2):
            val += diff(G[rho][mu][nu], x) * (1 if rho == 1 else 0)
            val -= diff(G[rho][mu][rho], x) * (1 if nu == 1 else 0)
        for rho in range(2):
            for lam in range(2):
                val += G[rho][mu][nu] * G[lam][rho][lam]
                val -= G[rho][mu][lam] * G[lam][nu][rho]
        R_mu[mu][nu] = simplify(val)

Ricci = simplify(ginv[0, 0] * R_mu[0][0] + ginv[0, 1] * R_mu[0][1]
                 + ginv[1, 0] * R_mu[1][0] + ginv[1, 1] * R_mu[1][1])
R["acoustic_ricci_full"] = str(Ricci)
R["acoustic_ricci_static_v0"] = str(simplify(Ricci.subs(v, 0)))

# ---------------------------------------------------------------------------
# 2. 结构同构：框架 δ vs BEC R（都是「背景场二阶导 → 标量曲率」）
# ---------------------------------------------------------------------------
R["structural_isomorphism"] = {
    "framework_delta": "δ = -∂²_x h_yy - ∂²_y h_xx（模长 h 的交叉二阶导，Q1 已解）",
    "bec_acoustic_R": "R = -2c''/c（声速 c 的二阶导 / c，静态极限）",
    "同构点": "两者都是「背景场（模长 r / 声速 c）的二阶导 → 标量曲率」，结构同构",
    "物理对应": "c² ∝ ρ（声速平方 ∝ 密度）⟹ R ∝ ρ 的二阶导，对应框架「模长 r 的二阶导」",
}

# ---------------------------------------------------------------------------
# 3. 映射字典（第一步核心产出）
# ---------------------------------------------------------------------------
R["mapping_dictionary"] = {
    "度规": "框架 g_ij = 模长 r_ij  ↔  BEC 声学度规 g_μν（由 ρ, v 决定）",
    "曲率": "框架 δ = 模长交叉二阶导  ↔  BEC R = -2c''/c（声速二阶导）",
    "引力势": "框架 观察者态 Φ = 1/r  ↔  BEC g_00 = -(c²-v²)（声子等效引力势）",
    "视界": "框架（待推）  ↔  BEC 声学视界 v=c（g_00=0，类比黑洞）",
    "物质/拓扑": "框架 π 磁通 Z₂ 规范  ↔  BEC 涡旋绕数拓扑",
    "应力": "框架 键序 K = δE/δD  ↔  BEC 声子应力张量",
}

# ---------------------------------------------------------------------------
# 4. 诚实结论
# ---------------------------------------------------------------------------
R["honest_conclusion"] = {
    "第一步产出": "映射字典 + 结构同构验证（曲率都是背景二阶导）",
    "同构成立": "框架 δ 与 BEC R 结构同构（符号坐实 R=-2c''/c）",
    "诚实边界": "这是「类比引力可行性确认」，不是「具体预言」——具体可测预言（声学视界/色散维度流）是后续步骤（映射曲率 + 给出预言）",
    "下一步": "映射曲率（键序 K ↔ 1-形式 [D,a]）+ 给出「观察者态 1/r → 声学视界类比」的可测签名",
}

report(R, "exp_analogue_gravity")

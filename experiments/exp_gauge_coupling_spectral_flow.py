"""
远景线 2（尺度读出）· 方案 C 第二刀（符号推导）：含谱流的完整 F² 项 = 0

第一刀错误（漏谱流）：只算 D_A^{-2} 的 V² × P_A 零阶，得 F² ∝ Tr[D^{-6}] ∝ 1/λ_min。
本脚本做完整推导（含谱流），符号验证「体项 + 边界项 = 0」，即 F² 项 = 0。

关键：f = (1/u)χ_{[λ_min, λ_c]}，F² 项 ∝ ∫ u f''(u) du。
f'' 含体项 + 边界 δ 项 + 边界 δ' 项。分部积分（分布意义）：
  ∫ u f''(u) du = [u f'(u) - f(u)]_0^∞ = 0（f 和 uf' 在边界 0,∞ 都 = 0）

本脚本验证：
  (1) 体项 ∫ (2/u³)χ · u du = 2(1/λ_min - 1/λ_c)；
  (2) 边界 δ 项 ∫ -(2/u²)[δ(u-λ_min)-δ(u-λ_c)] · u du = -2(1/λ_min - 1/λ_c)；
  (3) 体项 + 边界项 = 0（精确抵消）⟹ F² 项 = 0。
"""
from sympy import symbols, integrate, DiracDelta, simplify, Heaviside, oo, diff
from experiments._common import report

R = {}

u, lam_min, lam_c = symbols("u lambda_min lambda_c", positive=True)

# f(u) = (1/u)χ_{[λ_min, λ_c]} = (1/u)[θ(u-λ_min) - θ(u-λ_c)]
f = (1/u) * (Heaviside(u - lam_min) - Heaviside(u - lam_c))

# 体项：∫ (2/u³)χ · u du = ∫ (2/u²)χ du（在 [λ_min, λ_c] 内）
body = integrate(2/u**3 * u, (u, lam_min, lam_c))  # ∫ 2/u² du = 2(1/λ_min - 1/λ_c)
R["body_term"] = {
    "体项 ∫(2/u³)χ·u du": f"{simplify(body)} = 2(1/λ_min - 1/λ_c)",
}

# 边界 δ 项：∫ -(2/u²)[δ(u-λ_min)-δ(u-λ_c)]·u du = ∫ -(2/u)[δ(u-λ_min)-δ(u-λ_c)] du
# = -2(1/λ_min - 1/λ_c)（在边界取值）
boundary = integrate(-(2/u) * (DiracDelta(u - lam_min) - DiracDelta(u - lam_c)), (u, 0, oo))
R["boundary_term"] = {
    "边界 δ 项 ∫-(2/u)[δ(u-λ_min)-δ(u-λ_c)] du": f"{simplify(boundary)} = -2(1/λ_min - 1/λ_c)",
}

# 体项 + 边界项
R["cancellation"] = {
    "体项 + 边界项": f"{simplify(body + boundary)} = 0（精确抵消）",
    "F² 项": "= 0（硬截断谱作用量给不出规范场动力学）",
}

# 分部积分验证：∫ u f''(u) du = [u f'(u) - f(u)]_0^∞ = 0
# f 在 u=0（θ(0-λ_min)=0）和 u=∞（θ(∞-λ_min)-θ(∞-λ_c)=1-1=0）都 = 0
R["integration_by_parts"] = {
    "∫ u f''(u) du = [u f'(u) - f(u)]_0^∞": "f(0)=0（0<λ_min）、f(∞)=0（θ(∞-λ_min)=θ(∞-λ_c)=1）、uf' 边界也 = 0",
    "结论": "= 0，严格（分布意义的分部积分）",
}

R["honest_conclusion"] = {
    "第一刀错误": "只算 D_A^{-2} 的 V² × P_A 零阶，漏了谱流（谱投影 P_A 对 A 的依赖）。",
    "完整结果": "含谱流的 F² 项 = 0（体项和边界 δ 项精确抵消）——恢复了标准 CC 的 f4=f(0)=0 结论。",
    "物理含义": "硬截断 f=(1/u)χ 把低能谱（含规范场零模）整段截掉 ⟹ 规范场 F² 项消失，这是「结构结果」，不是「矩定义的问题」。",
    "方案 C 的结论": "从硬截断重新推导矩，不能救回规范场动力学——F² 项 = 0 是硬截断的结构边界。规范场动力学需别的来源（内涨落 1.11 B2、软截断、或新机制）。",
    "措辞": "完整符号推导（含谱流），推翻第一刀的 g²∝λ_min——F² 项 = 0，规范耦合发散是硬截断的结构结果。",
}

report(R, "exp_gauge_coupling_spectral_flow")

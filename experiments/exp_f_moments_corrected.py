"""
远景线 2（尺度读出）· f 来源定量修正：矩的物理对应（2026-09-27）

背景：exp_f_observer / exp_scale_readout / OQG 1.13 §9.2 的矩命名把物理对应标反了。
标准谱作用量（Chamseddine-Connes 1997 / Vassilevich）：
  S = Tr f(D^2/Lambda^2) = sum_k f_{2k} Lambda^{4-2k} a_{2k}
经 f(u) = int_0^oo e^{-su} phi(s) ds（Laplace 变换）：
  a0 阶（Lambda^4，宇宙学常数）: f0 = int f(u) u   du = int phi(s) s^{-2} ds
  a2 阶（Lambda^2，引力常数）:   f2 = int f(u)     du = int phi(s) s^{-1} ds
  a4 阶（Lambda^0，规范耦合）:   f4 = f(0)           = int phi(s) s^{0}  ds

错误（旧代码/OQG 1.13 §9.2）：
  「f0」= int f(u) u^{-1} du = 1/λ_min − 1/λ_c（红外幂次发散）—— 不是任何标准阶
  「f2」= int f(u) du       = ln(λ_c/λ_min)        —— 这才是标准 f2（引力）✓
  「f4」= int f(u) u du     = λ_c − λ_min           —— 这才是标准 f0（宇宙学常数），被错用在 g²

本脚本验证（修正版）：
  (1) 对照 f=e^{-u}：标准矩 f0=f2=f4=1（CC 1997）；
  (2) f = E = (1/u) chi_{[λ_min, λ_c]} 的正确矩；
  (3) 宇宙学常数 ρ_Λ = (λ_c − λ_min) Λ^4（O(1) 自然值）；
  (4) 规范耦合 f4 = f(0) = 0 的问题（硬截断伪影 → g² 发散）。
"""
from sympy import symbols, integrate, exp, oo, simplify, pi, ln
from experiments._common import report

R = {}

u = symbols("u", positive=True)

# ============ (1) 对照 f = e^{-u} ============
R["check_e_minus_u"] = {
    "f0 宇宙学常数 = int f u du": str(simplify(integrate(exp(-u) * u, (u, 0, oo)))),
    "f2 引力常数 = int f du": str(simplify(integrate(exp(-u), (u, 0, oo)))),
    "f4 规范耦合 = f(0)": str(simplify(exp(-u).subs(u, 0))),
    "意义": "CC 1997 标准：e^{-u} 给 f0=f2=f4=1，三个矩一致（坐实矩定义正确）",
}

# ============ (2) f = E = (1/u) chi 的正确矩 ============
lam_min, lam_c = symbols("lambda_min lambda_c", positive=True)
f0 = integrate((1/u) * u, (u, lam_min, lam_c))   # int f(u) u du = 宇宙学常数
f2 = integrate((1/u), (u, lam_min, lam_c))        # int f(u) du = 引力常数

R["correct_moments"] = {
    "宇宙学常数系数 f0 = int f(u) u du": f"{simplify(f0)} = λ_c − λ_min",
    "引力常数系数   f2 = int f(u) du": f"{simplify(f2)} = ln(λ_c/λ_min)",
    "规范耦合系数   f4 = f(0)": "0（硬截断：f 在 u<λ_min 处 = 0）",
    "对照旧错误 f0": "int f(u) u^{-1} du = 1/λ_min − 1/λ_c（红外幂次发散，非宇宙学常数）",
}

# ============ (3) 数值 N=128 ============
N_val = 128
lc_val = 2.0
lm_val = pi**2 / N_val**2
f0_num = lc_val - lm_val
f2_num = ln(lc_val / lm_val)

R["numerical_N128"] = {
    "λ_c = 2": lc_val,
    "λ_min = π²/N²": float(lm_val),
    "宇宙学常数系数 f0 = λ_c − λ_min": float(f0_num),
    "引力常数系数   f2 = ln(λ_c/λ_min)": float(f2_num),
    "宇宙学常数 ρ_Λ = f0 Λ^4 ≈ 2 Λ^4": "O(Λ^4) 自然值（层级 10^-122 未解决，但方向不再反）",
}

# ============ (4) 规范耦合 f(0) = 0 问题 ============
R["gauge_coupling_problem"] = {
    "f4 = f(0) = 0": "g² = 6π²/f4 = ∞（硬截断伪影）",
    "λ_min → 0 极限": "f(0) = lim_{u→0} 1/u = ∞ → g² = 0（另一个极端）",
    "结论": "f = (1/u)χ 的硬截断形式无法正确提取规范耦合——f(0) 对红外截断 λ_min 敏感（0 或 ∞），不是良定义。需更光滑的 f（软截断），或另取规范耦合的提取方式。",
    "旧结论影响": "「G 与规范耦合统一」（G/g² 纯数）在纠正后不成立——因为 f4=f(0)=0，G/g²=f4/(2π f2 Λ²)=0。该结论需撤回。",
}

R["honest_conclusion"] = {
    "发现": "旧代码/OQG 1.13 §9.2 把矩的物理对应标反了：宇宙学常数系数 = int f u du = λ_c−λ_min（O(1)），被错算成 int f u^{-1} du = 1/λ_min = N²/π²（红外发散、放大）。用户直觉「倒数错位」成立。",
    "纠正后": "宇宙学常数 ρ_Λ = (λ_c − λ_min) Λ^4 ≈ λ_c Λ^4 = 2 Λ^4（O(Λ^4) 自然值）；引力 G = 3π/(ln(λ_c/λ_min) Λ²)。",
    "新问题": "规范耦合 f4 = f(0) = 0（硬截断）→ g² 发散；「G 与规范耦合统一」结论需撤回。",
    "措辞": "修正矩对应（符号验证坐实），非「已解决宇宙学常数问题」——层级 10^-122 与规范耦合提取仍需新机制。",
}

report(R, "exp_f_moments_corrected")

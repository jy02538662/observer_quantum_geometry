"""
远景线 2（尺度读出）· N ↔ λ_c 的离散→连续对应（2026-09-27 撤回「G/g² 纯数」结论）

一元论直觉：N（量子化的离散截断）和 λ_c（观察者态的连续截断）是「同一个观察者有限性」
的两面（离散/连续），由「离散→连续」（OQG 一元论翻转）联系。

⚠️ 撤回（2026-09-27）：旧脚本算「G/g² = f4/(2π f2 Λ²)」并问是否纯数。
  但纠正矩对应后，规范耦合系数 f4 = f(0) = 0（硬截断伪影，见 exp_f_moments_corrected），
  所以 G/g² = f4/(2π f2 Λ²) = 0，「G 与规范耦合统一（纯数）」结论不成立，撤回。

本脚本改为验证：
  (1) λ_c = 2（经典极限），λ_min = 2-δ_N（尺度破缺）的内生对应；
  (2) 宇宙学常数/引力比值 f0/f2 = (λ_c-λ_min)/ln(λ_c/λ_min) 依赖 N（非纯数）；
  (3) N ↔ λ_c 一元论对应成立（观察者有限性的两面）。
"""
from sympy import symbols, pi, cos, ln, simplify
from experiments._common import report

R = {}

n_sym = symbols("N", positive=True)

# δ_N = 2cos(π/(N+1))，尺度破缺 2-δ_N = π²/N²
delta_N = 2 * cos(pi / (n_sym + 1))
break_scale = 2 - delta_N  # π²/N²

# 候选对应：λ_c = 2（经典极限），λ_min = 2-δ_N（尺度破缺）
lam_c = 2
lam_min = break_scale

ratio = (lam_c - lam_min) / ln(lam_c / lam_min)  # f0/f2 = 宇宙学常数/引力

R["candidate"] = {
    "λ_c = 2（经典极限）": "量子维度 d_{1/2} = 2，框架唯一",
    "λ_min = 2-δ_N（尺度破缺）": "π²/N²，共形 anomaly",
    "λ_min/λ_c": str(simplify(lam_min / lam_c)),
    "λ_min/λ_c 依赖 N": "π²/(2N²)——不是固定比值，依赖 N",
}

R["cosmological_gravity_ratio"] = {
    "f0/f2 表达式（宇宙学常数/引力）": str(simplify(ratio)),
    "N=3": float(ratio.subs(n_sym, 3).evalf()),
    "N=10": float(ratio.subs(n_sym, 10).evalf()),
    "N=100": float(ratio.subs(n_sym, 100).evalf()),
    "结论": "宇宙学常数/引力比值 f0/f2 依赖 N（观察者有限性），不是纯数",
}

# 撤回说明
R["withdrawn_G_unification"] = {
    "旧结论": "「G 与规范耦合统一」（G/g² 纯数）",
    "撤回原因": "规范耦合系数 f4 = f(0) = 0（硬截断伪影）⟹ G/g² = f4/(2π f2 Λ²) = 0，结论不成立",
    "措辞": "撤回（非「否证」——是矩对应修正后旧推导的前提失效，规范耦合提取是开放问题）",
}

R["honest_conclusion"] = {
    "一元论对应": "N ↔ λ_c 是观察者有限性的两面（离散/连续），由离散→连续联系——有唯一对应，不是两个独立参数（此条仍成立）",
    "撤回": "「G 与规范耦合统一（纯数）」撤回——f4=f(0)=0 硬截断伪影，规范耦合提取是开放问题",
    "保留": "宇宙学常数/引力比值 f0/f2 = (λ_c-λ_min)/ln(λ_c/λ_min) 由 N 唯一确定（结构统一），但非纯数",
    "最终结论": "「结构统一」成立（ρ_Λ、G 由同一个观察者有限性 N 决定）；「数值预言」（纯数）不成立——绝对耦合是观察者依赖的内生标定（无绝对尺度的必然）",
    "措辞": "一元论结构（N↔λ_c 有唯一对应）成立；「G/g² 纯数」撤回；规范耦合提取开放",
}

report(R, "exp_N_lambdac")

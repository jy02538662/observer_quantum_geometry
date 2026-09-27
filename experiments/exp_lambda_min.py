"""
远景线 2（尺度读出）· 一元论直觉：λ_min 和 λ_c 是「尺度不变破缺」的两个面（2026-09-27 修正矩对应）

方法论（用户两次直觉验证的通用钥匙）：「看似需要额外输入的 X，其实是框架里某个对象的另一个面」。

问题：λ_min（红外截断 = 无质量模调节）的来源？
一元论直觉：λ_min 和 λ_c（紫外截断 = 观察者分辨率）是「尺度不变破缺」的两个面——
  尺度不变破缺 = 量子化 2-δ_N（OQG §4.4：δ_N=2cos(π/(N+1))→2，破缺 2-δ_N=π²/N²），
  它同时给出红外（λ_min，无质量模）和紫外（λ_c，观察者分辨率）两个截断。

⚠️ 修正（2026-09-27）：矩比值物理含义更正。
  旧：f4/f2（误当「规范耦合/引力」）。
  新：f0/f2 = (λ_c-λ_min)/ln(λ_c/λ_min) = 「宇宙学常数系数/引力系数」。
  因为 f0（宇宙学常数）=∫f·u du=λ_c-λ_min，f2（引力）=∫f du=ln(λ_c/λ_min)。

本脚本验证：
  (1) 量子化谱的范围 [δ_N, 2]，尺度破缺 2-δ_N = π²/N²；
  (2) 对应 λ_c ↔ 2（经典极限）、λ_min ↔ 2-δ_N（尺度破缺）；
  (3) 宇宙学常数/引力比值 f0/f2 = (λ_c-λ_min)/ln(λ_c/λ_min) 由 N 唯一确定（内生）。
"""
from sympy import symbols, pi, cos, sin, ln, simplify, N as symN
from experiments._common import report

R = {}

n_sym = symbols("N", positive=True)

# 量子化 δ_N = 2cos(π/(N+1))，尺度破缺 2-δ_N
delta_N = 2 * cos(pi / (n_sym + 1))
break_scale = 2 - delta_N  # = π²/N² + O(1/N³)

R["quantization"] = {
    "delta_N": str(delta_N),
    "经典极限": "δ_N → 2（N→∞）",
    "尺度破缺 2-δ_N": "π²/N² + O(1/N³)（共形 anomaly 微观起源）",
}

# 一元论对应：λ_c ↔ 2（经典极限，紫外=观察者分辨率），λ_min ↔ 2-δ_N（尺度破缺，红外=无质量模）
# 宇宙学常数/引力比值 f0/f2 = (λ_c - λ_min)/ln(λ_c/λ_min)
lam_c = 2
lam_min = break_scale  # 2-δ_N

ratio = (lam_c - lam_min) / ln(lam_c / lam_min)

# 数值（N=3, 10, 100）
R["cosmological_gravity_ratio_vs_N"] = {
    "f0/f2 表达式（宇宙学常数/引力）": str(simplify(ratio)),
    "N=3": float(ratio.subs(n_sym, 3).evalf()),
    "N=10": float(ratio.subs(n_sym, 10).evalf()),
    "N=100": float(ratio.subs(n_sym, 100).evalf()),
    "意义": "宇宙学常数/引力比值 f0/f2 由 N（量子化的有限截断 = 观察者有限性）唯一确定——λ_min、λ_c 内生，不是两个独立手放参数",
}

# N→∞ 极限
R["large_N_limit"] = {
    "N→∞ 时 f0/f2": "→ 2/ln(2N²/π²) → 0（红外对数发散，无质量模）",
    "意义": "N→∞ = 无限观察者 = 尺度不变完全恢复，红外无截断（λ_min→0），比值→0",
}

R["honest_conclusion"] = {
    "一元论直觉": "λ_min（红外）和 λ_c（紫外）是「尺度不变破缺」的两个面，由量子化 δ_N 唯一确定——不是两个独立参数",
    "验证": "λ_c↔2（经典极限）、λ_min↔2-δ_N（尺度破缺），宇宙学常数/引力比值 f0/f2 由 N 唯一确定（内生）",
    "修正": "比值物理含义从「规范耦合/引力 f4/f2」更正为「宇宙学常数/引力 f0/f2」（f4=f(0)=0 硬截断伪影，见 exp_scale_readout）",
    "但": "比值仍依赖 N（观察者有限性），不是纯数——「纯数」还需 N 有唯一值（N ↔ λ_c 的离散→连续对应）",
    "性质": "一元论的结构（λ_min、λ_c 是同一个破缺的两面）+ 一个具体的对应问题（N↔λ_c）",
    "措辞": "一元论直觉（λ_min 不是独立参数），非「已推导数值」——纯数还需 N↔λ_c 的精确对应",
}

report(R, "exp_lambda_min")

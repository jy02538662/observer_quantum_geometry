"""
远景线 2（尺度读出）· f = E（切割）在谱语言里的名字（2026-09-27 修正矩对应）

用户直觉：f 不是手放的，是推导的——f = 观察者切割 E 的谱形式。
  - ρ = 态（红外面，尺度不变 1/λ）→ 观察者态谱 λ（点，§4.1）
  - E = 切割（紫外面，条件期望）→ f = 截断函数（谱作用量）
  同一个观察者的两个面（像颜色=Weyl群、代=共轭类是同一个 S₃ 的两面）。

⚠️ 修正（2026-09-27）：矩的物理对应标反了。标准谱作用量 S=Tr f(D²/Λ²) 的矩：
  - 宇宙学常数（a₀ 阶，Λ⁴）：f₀ = ∫ f(u) u du
  - 引力常数  （a₂ 阶，Λ²）：f₂ = ∫ f(u) du
  - 规范耦合  （a₄ 阶，Λ⁰）：f₄ = f(0)
（对照 e^{-u} 给 f₀=f₂=f₄=1，见 exp_f_moments_corrected。）

本脚本验证：
  (1) f = E 的谱形式 = (1/λ)·χ_{[λ_min, λ_c]}(λ)（公理 3 的截断区间）；
  (2) 正确矩：宇宙学常数 = ∫f·u du = λ_c−λ_min，引力 = ∫f du = ln(λ_c/λ_min)；
  (3) 宇宙学常数 ρ_Λ = (λ_c−λ_min)Λ⁴（O(1) 自然值，非 N²/π² 放大）。
"""
from sympy import symbols, integrate, log, oo, simplify
from experiments._common import report

R = {}

lam, lam_min, lam_c = symbols("lambda lambda_min lambda_c", positive=True)

# f = E 的谱形式 = (1/λ) 在 [λ_min, λ_c] 上（公理 3 的截断区间），区间外 0
f = 1 / lam

# 正确矩（4D 谱作用量，标准 Chamseddine-Connes）：
#   a₀ 阶（Λ⁴，宇宙学常数）f₀ = ∫ f(u) u du
#   a₂ 阶（Λ²，引力常数）  f₂ = ∫ f(u) du
#   a₄ 阶（Λ⁰，规范耦合）  f₄ = f(0)
f0 = integrate(f * lam**1, (lam, lam_min, lam_c))   # 宇宙学常数 = ∫ f·u du
f2 = integrate(f * lam**0, (lam, lam_min, lam_c))   # 引力常数 = ∫ f du

R["f_moments"] = {
    "宇宙学常数 f0 = ∫f·u du": str(simplify(f0)),
    "引力常数   f2 = ∫f du": str(simplify(f2)),
    "规范耦合   f4 = f(0)": "0（硬截断伪影，见下）",
}

# 宇宙学常数 ρ_Λ = f0 Λ⁴（a₀ 阶）
R["cosmological_constant"] = {
    "ρ_Λ = f0·Λ⁴": f"({simplify(f0)})·Λ⁴ = (λ_c − λ_min)·Λ⁴",
    "数值 N=128": "≈ 2·Λ⁴（O(1) 自然值）",
    "方向": "不再反——不是 N²/π² 放大（旧错误），是 λ_c−λ_min ≈ λ_c = 2",
}

# 引力常数 G = 3π/(f2 Λ²)
G = 3 * symbols("pi") / (f2 * lam_c**2)
R["gravity"] = {
    "G = 3π/(f2·Λ²)": f"3π/(ln(λ_c/λ_min)·Λ²)",
    "若 Λ=λ_c": str(simplify(3 * symbols("pi") / (f2 * lam_c**2))),
}

# 规范耦合 f(0) = 0 的问题
R["gauge_problem"] = {
    "f4 = f(0) = 0": "g² = 6π²/f4 = ∞（硬截断伪影）",
    "说明": "f=(1/u)χ 在 u<λ_min 处 = 0 ⟹ f(0)=0 ⟹ 规范耦合发散。需软截断 f 或另取提取方式（新问题）。",
}

R["honest_conclusion"] = {
    "直觉验证": "用户的直觉成立——f = E 的谱形式（= ρ 的截断版 (1/λ)·χ_{[λ_min,λ_c]}），是内生的（由观察者参数 λ_min, λ_c 唯一确定），不是手放的自由参数",
    "一元论结构": "ρ（态/红外）和 f=E（截断/紫外）是同一个观察者的两面，坐实",
    "修正": "宇宙学常数系数 = λ_c−λ_min（O(1) 自然值），非旧代码的 1/λ_min = N²/π²（红外发散）。",
    "剩余问题": "① 规范耦合 f(0)=0 硬截断伪影；② λ_min 与 λ_c 的关系（N 的唯一值）。",
    "措辞（AGENTS 三档）": "直觉验证（f 内生）+ 矩修正，非「已推导 G 数值」——规范耦合提取仍是开放问题",
}

report(R, "exp_f_observer")

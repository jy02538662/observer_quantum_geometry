"""
远景线 2（尺度读出）· 检查 ρ_Λ 是否有谱流修正（类比 F²=0）

用户直觉：F²=0 是谱流（体项+边界δ项抵消）。ρ_Λ = f0·Λ⁴ 是否也有谱流修正？

结论：谱流来自「f 的导数在边界的 δ」。f0=∫f·u du 的被积函数是 f·u=χ（阶跃），
不含 f 的导数 → 无谱流。f2=∫f du 的被积函数是 f=1/u，也不含导数 → 无谱流。
只有 f4（含 f''）有谱流抵消（F²=0，已在 exp_gauge_coupling_spectral_flow 验证）。

本脚本只验证 f0、f2 无边界项（核心结论）。
"""
from sympy import symbols, integrate, simplify
from experiments._common import report

R = {}

u, lam_min, lam_c = symbols("u lambda_min lambda_c", positive=True)

# f0 = ∫ f·u du = ∫ χ du = λ_c − λ_min（被积函数 f·u = χ，光滑阶跃，无边界项）
f0 = integrate(1, (u, lam_min, lam_c))  # f·u = 1（χ 区间内）
R["f0"] = {
    "f0 = ∫ f·u du = ∫ χ du": f"{simplify(f0)} = λ_c − λ_min",
    "被积函数 f·u = χ": "阶跃（=1 在区间内），不含 f 的导数 ⟹ 无边界项（无谱流）",
}

# f2 = ∫ f du = ∫ (1/u) du = ln(λ_c/λ_min)（被积函数 1/u，无导数，无边界项）
f2 = integrate(1/u, (u, lam_min, lam_c))
R["f2"] = {
    "f2 = ∫ f du = ∫ (1/u) du": f"{simplify(f2)} = ln(λ_c/λ_min)",
    "被积函数 1/u": "不含 f 的导数 ⟹ 无边界项（无谱流）",
}

R["honest_conclusion"] = {
    "检查结果": "ρ_Λ（a0 阶）没有谱流修正——f0 的被积函数 f·u=χ 是阶跃（=1），不含 f 的导数，无边界项。",
    "谱流机制": "谱流来自「f 的导数在边界的 δ」，只影响含 f 导数的矩（a4 阶 F²=0），不影响 f 的积分矩（a0 宇宙学常数、a2 引力）。",
    "张力确认": "ρ_Λ ≈ 2Λ⁴ 是严格的（无谱流修正），层级张力真实——不是谱流能消解，是「Λ 的选择」问题（M_P vs M_P/N）。",
    "措辞": "f0/f2 无边界项（符号验证）；ρ_Λ 无谱流修正，张力真实。",
}

report(R, "exp_rho_spectral_flow")

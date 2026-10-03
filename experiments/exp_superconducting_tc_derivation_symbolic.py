"""
符号验证：J_s = g·μ/(16π) 的完整推导链 + μ=Λ/λ_c 候选 A/B/C + μ² 坑

上一轮手推了相位刚度 J_s = ħ²n_s/(4m*) = g·μ/(16π)（2D Dirac 掺杂），
没验证。这里用 sympy 把每一步代数恒等式坐实，防止手推的因子错。

验证链：
  1. k_F = μ/(ħ v_F)                                  （Dirac 色散）
  2. n = g k_F²/(4π) = g μ²/(4π ħ² v_F²)              （载流子密度）
  3. m* = ħ k_F / v_F = μ / v_F²                      （Dirac 有效质量）
  4. J_s = ħ² n / (4 m*) = g μ / (16π)                （相位刚度，承重点①）
  5. T_BKT = (π/2) J_s = g μ / 32                     （BKT 转变温度）
  6. 候选 A/B/C：μ = Λ f(λ_c) 代入 T_BKT
  7. μ² 坑：填充 n ∝ μ² = 1/λ_c ⟹ μ = Λ/√λ_c（不是 Λ/λ_c）
"""
import sympy as sp
from experiments._common import report

R = {}

# ---- 符号 ----
mu, hbar, vF, g = sp.symbols("mu hbar vF g", positive=True)
Lambda, lam_c, lam_min = sp.symbols("Lambda lambda_c lambda_min", positive=True)

# ---- 1. Dirac 色散：k_F = μ/(ħ v_F) ----
kF = mu / (hbar * vF)

# ---- 2. 载流子密度（2D，g 重简并）----
n = g * kF**2 / (4 * sp.pi)
n_simplified = sp.simplify(n)

# ---- 3. Dirac 有效质量 ----
mstar = hbar * kF / vF
mstar_simplified = sp.simplify(mstar)

# ---- 4. 相位刚度 J_s = ħ² n / (4 m*) ----
Js = hbar**2 * n / (4 * mstar)
Js_simplified = sp.simplify(Js)

# ---- 5. BKT 转变温度 ----
T_BKT = (sp.pi / 2) * Js
T_BKT_simplified = sp.simplify(T_BKT)

R["chain"] = {
    "k_F = μ/(ħv_F)": str(kF),
    "n = g k_F²/(4π)": str(n_simplified),
    "m* = ħk_F/v_F": str(mstar_simplified),
    "J_s = ħ²n/(4m*)": str(Js_simplified),
    "T_BKT = (π/2)J_s": str(T_BKT_simplified),
}

# ---- 断言：J_s == g·μ/(16π)、T_BKT == g·μ/32 ----
assert sp.simplify(Js - g * mu / (16 * sp.pi)) == 0, "J_s ≠ gμ/16π"
assert sp.simplify(T_BKT - g * mu / 32) == 0, "T_BKT ≠ gμ/32"
R["assert_chain"] = {
    "J_s == g·μ/(16π)": "PASS",
    "T_BKT == g·μ/32": "PASS",
}

# ---- g=1 谷奇（笔记场景）与 g=2 自旋简并 ----
R["g_factor"] = {
    "g=1（谷奇 p+ip 无自旋，锁定自旋-谷）→ T_BKT = μ/32": str(sp.simplify(T_BKT.subs(g, 1))),
    "g=2（自旋简并 s 波）→ T_BKT = μ/16": str(sp.simplify(T_BKT.subs(g, 2))),
    "笔记用 μ/32 ⟺ g=1，必须显式写「谷奇锁定」这一步": True,
}

# ---- 6. 候选 A/B/C ----
mu_A = Lambda / lam_c
mu_B = Lambda * lam_min
mu_C = Lambda * (lam_c - lam_min)

Tc_A = sp.simplify(T_BKT.subs({g: 1, mu: mu_A}))
Tc_B = sp.simplify(T_BKT.subs({g: 1, mu: mu_B}))
Tc_C = sp.simplify(T_BKT.subs({g: 1, mu: mu_C}))

R["candidates"] = {
    "A: μ=Λ/λ_c → T_BKT": str(Tc_A),
    "B: μ=Λ·λ_min → T_BKT": str(Tc_B),
    "C: μ=Λ(λ_c−λ_min) → T_BKT": str(Tc_C),
}

# ---- 7. μ² 坑：若「填充比例 = 1/λ_c」且 n ∝ μ² ----
# n = g μ²/(4π ħ² v_F²) ∝ μ²；设 n/Λ² = 1/λ_c ⟹ μ = Λ/√λ_c
mu_sqrt = sp.symbols("mu_sqrt", positive=True)
# n ∝ μ²，n = const × μ²；「n ∝ 1/λ_c」⟹ μ² ∝ 1/λ_c ⟹ μ ∝ 1/√λ_c
R["mu2_pitfall"] = {
    "填充(粒子数) n ∝ μ²（Dirac 2D）": True,
    "若一元论对应是「填充=1/λ_c」⟹ μ² ∝ 1/λ_c ⟹ μ = Λ/√λ_c": "λ_c^{-1/2} 依赖",
    "若一元论对应是「能量 μ=偏好=1/λ_c」⟹ μ = Λ/λ_c": "λ_c^{-1} 依赖",
    "两个身份给不同函数形式（差根号），框架缺机制选定": True,
}

report(R, "exp_superconducting_tc_derivation_symbolic")

"""
继续推：定「掺杂 n ↔ 观察者 λ_c」是哪一对「两面」

上一轮结论：要锁死幂次，得回答「观察者有限性」用哪个量来量——
  「分辨率 1/λ_c」「模流能量 logλ_c」（已排除）「量子化数 N ∝ λ_min^{-1/2}」。

这一轮查框架自己的「观察者有限性」到底怎么定义的。框架里（f 来源 / 尺度读出）：
  f = E = (1/u)·χ_{[λ_min, λ_c]}(u)，u = D²/Λ²
  λ_min = π²/N²（量子化），λ_c = 2（UV 截断）
  f₀ = ∫f·u du = λ_c − λ_min（宇宙学常数矩）
  f₂ = ∫f du   = ln(λ_c/λ_min)（引力矩）

验证 + 对照：框架的「观察者有限性」是 f₀、f₂（谱矩），不是 1/λ_c（分辨率）。
这直接决定「掺杂 = 观察者有限性」该挂哪个量。
"""
import numpy as np
import sympy as sp
from experiments._common import report

R = {}

# ---- 1. 框架参数：λ_min = π²/N²，λ_c = 2 ----
N = 128
lam_min = np.pi**2 / N**2
lam_c = 2.0
R["framework_params"] = {
    "N（内部态数）": N,
    "λ_min = π²/N²": f"{lam_min:.6f}",
    "λ_c（UV 截断，框架值）": lam_c,
}

# ---- 2. 框架的「观察者有限性」度量：f₀、f₂ ----
f0 = lam_c - lam_min
f2 = np.log(lam_c / lam_min)
R["finiteness_moments"] = {
    "f₀ = λ_c − λ_min（宇宙学常数矩）": f"{f0:.6f}",
    "f₂ = ln(λ_c/λ_min)（引力矩）": f"{f2:.6f}",
    "1/λ_c（分辨率，笔记候选 A 用的）": f"{1/lam_c:.4f}",
    "关键": "f₀、f₂ 是「谱量」（λ_c→∞ 时 →∞），1/λ_c 才是「有限性」（λ_c→∞ 时 →0）",
}

# ---- 3. 符号验证：N ∝ λ_min^{-1/2}（量子化「数↔尺度」）----
lam_min_s, N_s = sp.symbols("lambda_min N", positive=True)
# λ_min = π²/N² ⟹ N = π/√λ_min ⟹ N ∝ λ_min^{-1/2}
N_expr = sp.pi / sp.sqrt(lam_min_s)
assert sp.simplify(N_expr - sp.pi * lam_min_s**(-sp.Rational(1, 2))) == 0
R["quantization_relation"] = {
    "λ_min = π²/N² ⟹ N = π/√λ_min": str(N_expr),
    "N ∝ λ_min^{-1/2}（数 ∝ 尺度^{-1/2}）": "符号 PASS",
    "注意": "这是 N ↔ λ_min（量子化面），不是 n ↔ λ_c（UV 面）",
}

# ---- 4. 对照：几种「观察者有限性」度量的 λ_c 依赖 ----
# 看「掺杂 = 观察者有限性」若挂不同量，给的 μ(λ_c) 幂次
def power_of(fn, lam_c_vals):
    """数值微分 d ln(fn) / d ln(λ_c)，看 λ_c 的幂次"""
    f = np.array([fn(lc) for lc in lam_c_vals])
    lf = np.log(f)
    llc = np.log(lam_c_vals)
    return np.polyfit(llc, lf, 1)[0]

lc_vals = np.array([1.5, 2.0, 2.5, 3.0, 4.0])
# 各度量（λ_min 固定 = 框架值）
measures = {
    "分辨率 1/λ_c": (lambda lc: 1.0 / lc, -1.0),
    "范围 f₀ = λ_c − λ_min": (lambda lc: lc - lam_min, None),
    "十年数 f₂ = ln(λ_c/λ_min)": (lambda lc: np.log(lc / lam_min), None),
}
R["power_compare"] = {}
for name, (fn, expected) in measures.items():
    p = power_of(fn, lc_vals)
    R["power_compare"][name] = {
        "d ln(度量)/d ln(λ_c)": round(float(p), 3),
        "含义": "λ_c 的幂次（-1=倒数，0=常数，+1=线性，对数≈0 但缓慢）",
    }

# ---- 5. 净结论 ----
R["conclusion"] = {
    "纠正（上一版结论写反了）": "f₀=λ_c−λ_min、f₂=ln(λ_c/λ_min) 是「谱量」（λ_c→∞ 时 →∞），不是「有限性」。",
    "「观察者有限性」正确定义": "λ_c→∞ = 全知 = 无限性 → 0 的量：1/λ_c（分辨率）、λ_min/λ_c（比值）——这些才是「有限性」。",
    "笔记候选 A 用 1/λ_c（分辨率）": "方向对（λ_c→∞ 时 →0 = 无限性）。",
    "量子化「数↔尺度」N∝λ_min^{-1/2}": "是 N↔λ_min（IR/量子化面），不是 n↔λ_c（UV 面）——之前当作「掺杂↔λ_c」的候选 λ_c^{-1/2} 是红鲱鱼。",
    "所以「两面」对": "「掺杂 ↔ 观察者有限性」= 「填充 n ↔ 分辨率 1/λ_c」这一对，不是「数↔尺度」（那是 N↔λ_min）。",
    "剩的最后一步": "「n ∝ 1/λ_c」本身仍是 ansatz（两面给了框架：同一个量两个投影，但没说具体变换「n = 1/λ_c」）。候选 c（μ=Λ/√λ_c）是剩下的候选，但精确的「有限性」量（1/λ_c vs λ_min/λ_c vs 1/ln）还没锁死。",
}

report(R, "exp_superconducting_doping_finiteness")

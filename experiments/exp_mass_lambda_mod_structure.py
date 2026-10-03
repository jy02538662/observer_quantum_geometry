"""
质量谱动力学 · λ_mod 与 λ_c 的关系（独立推导：从 ρ=C/λ 推模流频率）

用户澄清：λ_mod 不是截断，是「模流频率」；截断是 λ_c（UV）和 λ_min（IR）。
要推的关系：λ_mod 和 λ_c 到底怎么连。

独立推导（不套笔记公式，从观察者态 ρ=C/λ 直接算模流频率）：
  1. 观察者态 ρ = C/λ（λ∈[λ_min,λ_c]），归一化 C = 1/ln(λ_c/λ_min)。
  2. 模流生成元 = log ρ = log C − log λ。
  3. 模流频率 λ_mod = log ρ 谱的最大特征值 = log C − log λ_min（λ_min 处最大）。
  4. 展开：λ_mod = log(1/λ_min) − log(ln(λ_c/λ_min))。

关键结构：λ_mod 主要由「IR 截断 λ_min 的倒数对数」决定，λ_c 只通过「对数宽度
ln(λ_c/λ_min)」进一个次领头的「对数的对数」。

数值（N=128，λ_min=π²/N²，λ_c=2）：
  log(1/λ_min) = log(N²/π²)
  ln(λ_c/λ_min) = ln(2N²/π²)
  λ_mod = log(N²/π²) − log(ln(2N²/π²)) = log(N²/(π² ln(2N²/π²)))  ← 与 exp_mass_final 一致
"""
import numpy as np
import sympy as sp
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. 观察者态 ρ = C/λ（归一化）
# ---------------------------------------------------------------------------
lam, lam_c, lam_min = sp.symbols("lambda lambda_c lambda_min", positive=True)
# 归一化：∫_{λ_min}^{λ_c} ρ dλ = 1 → C ∫ dλ/λ = C ln(λ_c/λ_min) = 1 → C = 1/ln(λ_c/λ_min)
C = 1 / sp.log(lam_c / lam_min)
rho = C / lam

R["observer_state"] = {
    "观察者态 ρ = C/λ": str(rho),
    "归一化 C = 1/ln(λ_c/λ_min)": str(C),
    "λ ∈ [λ_min, λ_c]": "有限区间（IR 截断 λ_min，UV 截断 λ_c）",
}

# ---------------------------------------------------------------------------
# 2. 模流生成元 log ρ，最大特征值 = λ_mod
# ---------------------------------------------------------------------------
log_rho = sp.log(C / lam)  # = log C − log λ
# 最大特征值：λ 最小（λ_min）处，log ρ 最大
lam_mod_expr = sp.log(C / lam_min)

R["modular_flow"] = {
    "模流生成元 log ρ = log C − log λ": str(log_rho),
    "模流频率 λ_mod = log(C/λ_min)（λ_min 处最大）": str(lam_mod_expr),
}

# ---------------------------------------------------------------------------
# 3. 展开 λ_mod = log(1/λ_min) − log(ln(λ_c/λ_min))
# ---------------------------------------------------------------------------
# log(C/λ_min) = log(1/(λ_min ln(λ_c/λ_min))) = log(1/λ_min) − log(ln(λ_c/λ_min))
term_IR = sp.log(1 / lam_min)          # IR 截断的倒数对数
term_width = sp.log(sp.log(lam_c / lam_min))  # 对数宽度的对数
expanded = term_IR - term_width

R["expansion"] = {
    "λ_mod = log(1/λ_min) − log(ln(λ_c/λ_min))": str(expanded),
    "IR 项 log(1/λ_min)": "主导（IR 截断的倒数对数）",
    "宽度项 log(ln(λ_c/λ_min))": "次领头（λ_c 只进这里，且是对数的对数）",
    "符号验证展开 = 原式": str(sp.simplify(expanded - lam_mod_expr) == 0),
}

# ---------------------------------------------------------------------------
# 4. 数值（N=128）
# ---------------------------------------------------------------------------
def lam_mod_num(N, lam_c_val=2.0):
    lam_min_val = np.pi**2 / N**2
    term_ir = np.log(1 / lam_min_val)
    term_w = np.log(np.log(lam_c_val / lam_min_val))
    return lam_min_val, term_ir, term_w, term_ir - term_w

for N in [64, 128, 256]:
    lm, t_ir, t_w, lm_mod = lam_mod_num(N)
    R[f"N={N}"] = {
        "λ_min = π²/N²": f"{lm:.3e}",
        "IR 项 log(1/λ_min)": round(float(t_ir), 3),
        "宽度项 log(ln(λ_c/λ_min))": round(float(t_w), 3),
        "λ_mod = IR 项 − 宽度项": round(float(lm_mod), 4),
    }

# 与 exp_mass_final 的 λ_mod 公式对齐
N = 128
lm_final = np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))
R["align_with_mass_final"] = {
    "λ_mod（独立推导）= log(N²/π²) − log(ln(2N²/π²))": round(float(np.log(N**2/np.pi**2) - np.log(np.log(2*N**2/np.pi**2))), 4),
    "λ_mod（exp_mass_final）= log(N²/(π² ln(2N²/π²)))": round(float(lm_final), 4),
    "两者一致": bool(np.isclose(np.log(N**2/np.pi**2) - np.log(np.log(2*N**2/np.pi**2)), lm_final, atol=1e-9)),
}

R["honest_conclusion"] = {
    "独立推导结果": "λ_mod = log(C/λ_min) = log(1/λ_min) − log(ln(λ_c/λ_min))——模流频率 = 「IR 截断的倒数对数」减「对数宽度的对数」。",
    "关键结构": "λ_mod 由 IR 截断 λ_min 主导（log(1/λ_min)），λ_c 只通过「对数的对数」进次领头项。所以 λ_mod 和 λ_c **不是「两面」**，是「IR 截断的对数放大，λ_c 只做参考宽度」。",
    "⚠️ 笔误修正": "log(1/λ_min) = log(N²/π²) ≈ 7.41（不是 8.11）；8.11 = ln(λ_c/λ_min) = ln(2N²/π²)。结构不变。",
    "对质量指数的意义": "m_n = e^{λ_mod·n} 的「指数 λ_mod」=「IR 截断（尺度破缺 2−δ_N）的对数放大」——即「有限性（IR 截断）→ 对数 → 指数层级」。这坐实了「质量指数来自穿过有限性（隧穿）」：IR 截断 = 有限性，对数放大 = 隧穿指数。",
    "下一步": "把「λ_mod = IR 截断的对数」接到「绕数 = 代阶」（n 接拓扑荷），完成「有限性 → 指数 → 拓扑荷」闭环。",
    "定位": "【独立推导，符号坐实】。λ_mod 与 λ_c 的关系精确了：不是两面，是「IR 主导 + UV 次领头」。",
}

report(R, "exp_mass_lambda_mod_structure")

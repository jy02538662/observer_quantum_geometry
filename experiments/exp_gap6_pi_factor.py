"""
v11 缺口 6 · c 的精确值：π 因子落在哪

问题：ρ_Λ = 2(M_P/N_ext)^4 匹配 HDE ρ_Λ = 3c² M_P²/R²，
c = √(2/3)（无 π，w=-1.016，DESI 一致）vs √(2/(3π))（有 π，w=-1.54，矛盾）。

精确化链条，定位 π 因子：
  1. λ_min = π²/N²（量子化尺度破缺，精确：2 - 2cos(π/(N+1)) ≈ π²/(N+1)²）
  2. λ_min = L_P/R = 1/(R M_P)（识别：红外截断 = 普朗克尺度/宇宙尺度，无 π）
  3. ⟹ N = π √(R M_P)（精确推导，带 π）
  4. ρ_Λ = 2 M_P^4 / N^4
  5. HDE 匹配 ⟹ c²

三种 N 定义给三个 c：
  A. N = √(R M_P)         （无 π，量级对应）        → c² = 2/3
  B. N = π^{1/4} √(R M_P) （Bekenstein-Hawking 熵） → c² = 2/(3π)
  C. N = π √(R M_P)       （量子化精确）            → c² = 2/(3π^4)

结论：π 因子落在「量子化弦数 N（带 π²）→ 宇宙学熵 S=N^4（无 π）」的映射上。
观测（DESI 倾向 w=-1.016）选 A（无 π），但这是「观测选边」非「框架逼出」。
"""
import numpy as np
from experiments._common import report

R = {}

# 量子化尺度破缺 λ_min 的精确展开
# λ_min = 2 - δ_N = 2(1 - cos(π/(N+1))) ≈ π²/(N+1)²
import sympy as sp
N_sym = sp.symbols('N', positive=True)
lam_min_exact = 2 - 2*sp.cos(sp.pi/(N_sym+1))
lam_min_expand = sp.series(lam_min_exact, N_sym, sp.oo, 5)
R["lambda_min"] = {
    "精确 λ_min = 2(1-cos(π/(N+1)))": str(lam_min_exact),
    "大 N 展开": str(lam_min_expand),
    "领头阶": "π²/(N+1)² ≈ π²/N²（带 π²）",
}

# 三个 N 定义 → 三个 c²
c2_A = 2/3                    # N = √(RM_P)，无 π
c2_B = 2/(3*np.pi)            # N = π^{1/4}√(RM_P)，BH 熵
c2_C = 2/(3*np.pi**4)         # N = π√(RM_P)，量子化精确

def w_of_c(c):
    Omega0 = 0.7
    return -(1/3) * (1 + 2*np.sqrt(Omega0)/c)

R["three_c"] = {
    "A 无π（N=√RM_P）": f"c=√(2/3)={np.sqrt(c2_A):.4f}, w={w_of_c(np.sqrt(c2_A)):.4f}",
    "B BH熵（N=π^1/4√RM_P）": f"c=√(2/(3π))={np.sqrt(c2_B):.4f}, w={w_of_c(np.sqrt(c2_B)):.4f}",
    "C 量子化精确（N=π√RM_P）": f"c=√(2/(3π^4))={np.sqrt(c2_C):.4f}, w={w_of_c(np.sqrt(c2_C)):.4f}",
    "DESI 2024": "w0 ≈ -0.99 ± 0.06，倾向 A（无 π）",
}

# π 因子落在哪个环节：量化
# λ_min = π²/N² = 1/(R M_P) ⟹ N = π√(R M_P)
# 而 ρ_Λ = 2M_P^4/N^4，用 N=π√(RM_P) 给 c²=2/(3π^4)
# 若用 N=√(RM_P)（丢 π）给 c²=2/3
R["pi_location"] = {
    "量子化 λ_min 的 π²": "λ_min = π²/N²（Chebyshev 零点的 π，真实）",
    "宇宙学 λ_min = L_P/R": "无 π（L_P=1/M_P 自然单位，R 视界半径）",
    "识别 ⟹ N = π√(RM_P)": "带 π —— 弦数 N 与宇宙学 N_ext=√(RM_P) 差 π 倍",
    "π 的物理候选": "一维弦 → 二维面积 → 四维熵的维度映射（圆周率），未证",
}

R["honest_conclusion"] = {
    "精确推": "若严格走 λ_min=π²/N² = L_P/R，得 N=π√(RM_P)，c²=2/(3π^4)，c≈0.083，w≈-7 —— 荒谬，说明链条是「量级对应」非「精确等式」。",
    "c=√(2/3) 的前提": "S=R²M_P²（无 π）的量级定义，即 N_ext=√(RM_P)（丢量子化的 π）。",
    "缺口 6 的诚实答案": "π 因子是「量子化弦数（带 π²）→ 宇宙学熵（无 π）」映射的待定系数，无法从第一性原理逼出。观测（DESI）选边「无 π」（c=√(2/3)），框架暂给不出无 π 的第一性理由。这是「待定系数 + 观测选边」，不是「结构矛盾」也不是「算错」。",
    "措辞": "三档「未证实」——π 因子是待定系数（默认档），非「否证」（没证明 π 因子绝不可能有框架理由），也非「重言式」。",
}

report(R, "exp_gap6_pi_factor")

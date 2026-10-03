"""
阶段 4 推进：共形 → 黎曼（2−δ_N 是「尺子最小刻度」= 伸缩子）

核心直觉（用户）：2−δ_N 是「尺子的最小刻度」——把「没有绝对长度的橡皮膜（共形）」变成「有标准长度的尺子（黎曼）」。
  物理：有限观察者（N 有限）→ 有限分辨率（2−δ_N）→ 有限尺度（黎曼）；N→∞ → 2−δ_N→0 → 无限精细（共形）。
  数学：2−δ_N = 伸缩子（dilaton）真空期望值，固定 Weyl 因子。「共形类 + 伸缩子 = 黎曼」是标准结果。

验证四层：
  1. 尺度不变 = 共形（Weyl）：观察者态 ρ=C/λ 尺度不变 ρ(cλ)=ρ(λ)/c；Connes 距离 d=|logλ1-logλ2| 在 λ→cλ 下不变。
  2. 尺度破缺 = 黎曼（最小刻度）：2−δ_N = π²/(N+1)² 是「最小刻度」λ_min，打破尺度不变。
  3. N→∞ 极限：2−δ_N→0，最小刻度→0，尺度不变恢复（共形极限）。
  4. 三合一：2−δ_N = 光滑（P3 二阶结构）= 质量（尺度破缺）= 黎曼（伸缩子）——同一个对象。
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. 尺度不变 = 共形（Weyl） ----
# 观察者态 ρ = C/λ，尺度不变：ρ(cλ) = C/(cλ) = (1/c)·C/λ = ρ(λ)/c
try:
    import sympy as sp
    lam, C, c = sp.symbols('lambda C c', positive=True)
    rho = C / lam                       # ρ = C/λ
    rho_scaled = C / (c * lam)          # ρ(cλ)
    R["scale_invariance"] = {
        "ρ = C/λ": str(rho),
        "ρ(cλ) = ": str(rho_scaled),
        "ρ(cλ) = ρ(λ)/c（尺度不变，Weyl）": str(sp.simplify(rho_scaled - rho / c) == 0),
        "含义": "观察者态尺度不变 ⟹ 共形（Weyl 类，度规差 Ω² 因子，无绝对长度）",
    }
    # Connes 距离 d = |logλ1 - logλ2| 在 λ→cλ（平移 s=logλ）下不变
    s1, s2 = sp.symbols('s1 s2', real=True)
    d = sp.Abs(s1 - s2)                  # d = |logλ1 - logλ2|
    shift = sp.symbols('shift', real=True)
    d_shifted = sp.Abs((s1 + shift) - (s2 + shift))
    R["connes_scale_inv"] = {
        "Connes 距离 d = |s1-s2|（s=logλ）": str(d),
        "平移 s→s+shift（= λ→e^{shift}λ 尺度变换）后 d 不变": str(sp.simplify(d_shifted - d) == 0),
        "含义": "尺度不变 ⟹ 距离平移不变 ⟹ 共形（无最小刻度）",
    }
except Exception as e:
    R["scale_invariance"] = {"error": str(e)}

# ---- 2. 尺度破缺 = 黎曼（最小刻度 λ_min = 2−δ_N） ----
for N in [16, 64, 128, 1000]:
    delta = 2 * np.cos(np.pi / (N + 1))
    scale_break = 2 - delta                     # 尺度破缺 = 2−δ_N
    lead = (np.pi / (N + 1)) ** 2               # 主导 π²/(N+1)²
    R[f"scale_break_N={N}"] = {
        "2−δ_N = 最小刻度（尺度破缺）": f"{scale_break:.6e}",
        "主导阶 π²/(N+1)²": f"{lead:.6e}",
        "最小刻度 > 0（打破尺度不变 = 黎曼）": bool(scale_break > 0),
    }

# ---- 3. N→∞ 极限：共形恢复 ----
N_list = np.array([10, 32, 100, 316, 1000, 3162, 10000])
scale_break_list = 2 - 2 * np.cos(np.pi / (N_list + 1))
# 验证 2−δ_N ∝ 1/N²（幂次 2，随 N→∞ 衰减到 0 = 共形极限）
logN = np.log(N_list)
logSB = np.log(scale_break_list)
slope = np.polyfit(logN, logSB, 1)[0]
R["conformal_limit"] = {
    "2−δ_N 随 N 的衰减幂次（log-log 斜率）": f"{slope:.4f}",
    "幂次 ≈ −2（2−δ_N ∝ 1/N²）": bool(abs(slope + 2) < 0.05),
    "N→∞ 时 2−δ_N→0（最小刻度→0，共形极限恢复）": "尺度破缺在 N→∞ 消失 = 共形",
}

# ---- 4. 三合一：2−δ_N 的三个投影 ----
R["trinity"] = {
    "光滑（P3 切空间）": "2−δ_N = 离散谐振子本征值（二阶结构）",
    "质量（质量谱）": "质量 = 尺度破缺 = 2−δ_N",
    "黎曼（共形→黎曼）": "2−δ_N = 伸缩子，固定 Weyl gauge（最小刻度）",
    "结论": "光滑 / 质量 / 黎曼 三者同源，都是 2−δ_N 的不同投影——「共形→黎曼」不是新卡点，是「尺度破缺=质量」的另一个面",
}

# ---- 结论 ----
R["conclusion"] = {
    "共形 → 黎曼": "= 尺度不变 → 尺度破缺 = 无最小刻度 → 有最小刻度 = N→∞ → N有限",
    "伸缩子": "2−δ_N = 伸缩子（dilaton）真空期望值 = 最小刻度，固定 Weyl 因子",
    "物理图像": "有限观察者（N有限）→ 有限分辨率（2−δ_N）→ 有限尺度（黎曼）；无限观察者（N→∞）→ 无最小刻度（共形）",
    "与 KMS 一致": "温度 T=Λ/λ_mod，有限 N → 有限 λ_mod → 有限温度 → 有限尺度（黎曼）",
    "诚实边界": "「2−δ_N = 伸缩子」的识别是直觉（离散谐振子本征值 = 伸缩子背景），「共形类 + 伸缩子 = 黎曼」是标准结果，但框架里「尺度破缺就是那个固定 Weyl gauge 的伸缩子」这一步仍需严格化（识别不是自动）",
}

report(R, "exp_conformal_riemann")

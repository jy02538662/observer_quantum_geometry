"""
继续推：从谱维数 d_s=2 导出态密度 ν(E) ∝ E（框架内推导 n ∝ μ²）

上一轮「载流子密度 n ∝ μ²」是手数态（格点数态）验证的，没有回答「为什么 ∝ μ²」。
这一轮坐实「为什么」：**谱维数 d_s=2 ⟹ ν(E)∝E ⟹ n∝μ²**。

链（标准热核 → 态密度关系）：
  d_s = 2 ⟺ Tr(e^{-tD²}) ∝ t^{-d_s/2} = 1/t  (t→∞, 低能)
  Tr(e^{-tD²}) = ∫ν(E) e^{-tE²} dE  ⟹  1/t ⟺ ν(E) ∝ E  (逆 Laplace)
  ν(E) ∝ E ⟹ n(μ) = ∫_0^μ ν(E)dE ∝ μ²

格点 Dirac 模型 H = sin k_x σ_x + sin k_y σ_y（4 个 Dirac 点，谱 ε_±=±√(sin²k_x+sin²k_y)）。
验证：
  A. 热核迹 Tr(e^{-tD²}) ∝ 1/t（谱维数 d_s=2）
  B. 态密度 ν(E) ∝ E（线性，低能）
  C. n(μ) = ∫ν dE ∝ μ²（= 上一轮载流子密度的「为什么」）
"""
import numpy as np
from experiments._common import report

R = {}


def dirac_bands(L):
    ks = (2 * np.pi / L) * np.arange(L) - np.pi
    KX, KY = np.meshgrid(ks, ks)
    e2 = np.sin(KX) ** 2 + np.sin(KY) ** 2   # ε_±² = sin²kx + sin²ky
    return e2.flatten()


# ---- A. 热核迹 Tr(e^{-tD²}) ∝ 1/t ----
L = 600
e2 = dirac_bands(L)  # ε_±² 的平方（两带共用，故乘 2）

ts = np.array([10.0, 20.0, 40.0, 80.0, 160.0])
heat = np.array([np.sum(2 * np.exp(-t * e2)) / (L * L) for t in ts])  # 2 带 × 每带 e^{-tε²}

# 对数斜率：log(heat) vs log(t)，斜率应 = -1（d_s=2）
log_t = np.log(ts)
log_h = np.log(heat)
slope = np.polyfit(log_t, log_h, 1)[0]
R["A_heat_kernel"] = {
    "t 扫描": {f"t={int(t)}": round(float(heat[i]), 6) for i, t in enumerate(ts)},
    "对数斜率 d log(Tr)/d log(t)": round(float(slope), 4),
    "理论 −d_s/2 = −1（d_s=2）": -1.0,
    "谱维数 d_s=2（热核 ∝ 1/t）": bool(abs(slope + 1.0) < 0.05),
}

# ---- B. 态密度 ν(E) ∝ E（直方图，低能线性） ----
L2 = 600
ks2 = (2 * np.pi / L2) * np.arange(L2) - np.pi
KX2, KY2 = np.meshgrid(ks2, ks2)
eps_pos = np.sqrt(np.sin(KX2) ** 2 + np.sin(KY2) ** 2)  # 上带 ε_+ = |E|
E_vals = eps_pos.flatten()

bins = 120
E_max = 0.4
hist, edges = np.histogram(E_vals, bins=bins, range=(0, E_max))
centers = (edges[:-1] + edges[1:]) / 2

# 低能段线性拟合：ν(E) ≈ a·E，检验 a 非零且线性
mask = centers < 0.25
lin = np.polyfit(centers[mask], hist[mask], 1)
R["B_dos_linear"] = {
    "低能段拟合 ν(E) ≈ a·E + b": f"a={lin[0]:.1f}, b={lin[1]:.1f}",
    "理论 ν(E) ∝ E（b≈0，线性）": "a>0 且 b 相对小",
    "DOS 线性（2D Dirac）": bool(lin[0] > 0 and abs(lin[1]) < 5),
}

# ---- C. n(μ) = ∫ν dE ∝ μ²（从 DOS 积分，非手数态） ----
# 累积：n(μ) = Σ_{E<μ} ν(E)（连续 = 累积直方图）
mus = np.array([0.1, 0.2, 0.3, 0.4])
n_cum = np.array([np.sum(E_vals < m) / (L2 * L2) for m in mus])
theory_n = mus**2 / np.pi   # 4 个 Dirac 点 × μ²/(4π)

R["C_dos_integral"] = {
    "μ 扫描 n(μ)=Σ_{E<μ}": {f"{m}": round(float(n_cum[i]), 6) for i, m in enumerate(mus)},
    "理论 μ²/π": {f"{m}": round(float(theory_n[i]), 6) for i, m in enumerate(mus)},
    "n ∝ μ²（从 DOS 积分，= 谱维数 d_s=2 的推论）": bool(np.max(np.abs(n_cum - theory_n) / theory_n) < 0.03),
}

# ---- 净结论 ----
R["conclusion"] = {
    "坐实的": "① 谱维数 d_s=2（热核 ∝ 1/t，斜率 -1）；② 态密度 ν(E) ∝ E（低能线性）；③ n(μ) ∝ μ²（从 ν 积分）。三者闭环 = 「n∝μ²」的框架内推导（不是手数态）。",
    "剩最后一步": "n ∝ 1/λ_c（「填充 = 观察者有限性」的对应 + 幂次）仍是 ansatz——这是「缺态密度 ν 的 λ_c 依赖」的物质侧墙，见笔记 [[掺杂μ↔观察者截断λ_c：推导尝试（身份=粒子数，缺态密度ν）]]。",
    "若 n ∝ 1/λ_c 成立": "μ = Λ/√λ_c（候选 c，不是笔记的 Λ/λ_c）。",
}

report(R, "exp_superconducting_dos_spectral_dim")

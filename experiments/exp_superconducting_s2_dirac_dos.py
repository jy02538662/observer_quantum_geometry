"""
验证：S² 上的 Dirac 算子谱 + 态密度 ν₂(E) ∝ E（「秤」在哪）

思路（用户提出）：框架的「内部 SU(2) → S²」自带一个 2D Dirac 算子，
它的态密度 ν₂(E) ∝ E 就是「数粒子」的秤（不是用 ρ=C/λ 去数）。

S²（半径 R=1）上自旋 1/2 Dirac 算子的谱是标准解析结果：
  本征值 E = ±(k+1)，k = 0,1,2,...，每个符号的简并度 2(k+1)
  （即 ±1 ×2、±2 ×4、±3 ×6、±4 ×8、...）

验证：
  A. 构造谱（解析），算 DOS，验证 ν₂(E) ∝ E（线性）
  B. n(μ) = ∫_0^μ ν₂ dE ∝ μ²（= 「掺杂粒子数」的幂律，非对数）
  C. 诚实标注：这把秤给的是「n ∝ μ²」，但「n ∝ 1/λ_c」仍不是这把秤给的
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 构造 S² Dirac 谱（解析）----
kmax = 2000
eigenvals = []
for k in range(kmax):
    val = k + 1                      # 正能级 E = k+1
    deg = 2 * (k + 1)                # 简并度 2(k+1)
    eigenvals.extend([val] * deg)    # 正能谱
    eigenvals.extend([-val] * deg)   # 负能谱（粒子-空穴对称）
E = np.array(eigenvals)

# ---- A. DOS ν₂(E) ∝ E ----
Emax = 200
bins = 80
hist, edges = np.histogram(E, bins=bins, range=(0, Emax))
centers = (edges[:-1] + edges[1:]) / 2
lin = np.polyfit(centers[centers < 150], hist[centers < 150], 1)
R["A_dos_linear"] = {
    "谱：±(k+1)，简并 2(k+1)（S² 标准结果）": True,
    "DOS 拟合 ν₂(E) ≈ a·E + b": f"a={lin[0]:.2f}, b={lin[1]:.2f}",
    "ν₂(E) ∝ E（线性，b≈0）": bool(lin[0] > 0 and abs(lin[1]) < 50),
}

# ---- B. n(μ) = ∫ν dE ∝ μ² ----
# 累积态数（正能，能量 < μ）：N(μ) = Σ_k 2(k+1)，k 从 0 到 μ-1
mus = [20.0, 40.0, 60.0, 80.0]
N_cum = []
for mu in mus:
    k_occ = int(np.floor(mu)) - 1
    N = 2 * (k_occ + 1) * (k_occ + 2) / 2   # Σ_{k=0}^{k_occ} 2(k+1) = (k_occ+1)(k_occ+2)
    N_cum.append(N)
N_cum = np.array(N_cum)
mus_np = np.array(mus)
# 拟合 N ∝ μ^p，p 应 ≈ 2
log_mu = np.log(mus_np)
log_N = np.log(N_cum)
p = np.polyfit(log_mu, log_N, 1)[0]
R["B_n_quadratic"] = {
    "μ 扫描 N(μ)": {f"{int(m)}": round(float(N_cum[i]), 1) for i, m in enumerate(mus)},
    "拟合 N ∝ μ^p 的 p": round(float(p), 3),
    "理论 p=2（n ∝ μ²）": 2.0,
    "n ∝ μ²（S² Dirac 直接给）": bool(abs(p - 2.0) < 0.05),
}

# ---- C. 诚实标注：这把秤给什么、不给什么 ----
R["C_honest"] = {
    "秤给的东西": "S² 上的 Dirac 态密度 ν₂(E) ∝ E ⟹ n(μ) ∝ μ²（幂律，非对数）。这是「数粒子」的正确秤，不是 ρ=C/λ。",
    "秤没给的东西": "「n ∝ 1/λ_c」——把「填充」接到「观察者有限性 λ_c」这一步，S² Dirac 没给。这把秤数「每个能量多少个态」，但没说「填充多少 = 观察者多有限」。",
    "净结论": "S² Dirac 把「n ∝ μ²」从「手数态」升级成「内部 SU(2)→S² 结构自带」，这是真进展（纠正「框架没数对象」）。但 μ = Λ/√λ_c 仍差「n ∝ 1/λ_c」这最后一步，它不是 S² Dirac 能给的。",
}

report(R, "exp_superconducting_s2_dirac_dos")

"""
离散圈图 · 谱投影积掉高能模，看有效耦合 g²(λ′) 的 λ′ 依赖（跑动）

方法（Wilsonian 谱投影，不撞主墙）：
  1. 构造 π 磁通 D(L)，加局域通量 δφ；
  2. 对角化得本征值 λ_k(δφ)；
  3. 硬截断谱投影：只保留 λ_k < λ′（积掉高能模 λ_k > λ′）；
  4. 有效谱和 S_eff(δφ, λ′) = Σ_{λ_k < λ′} e^{-λ_k²/Λ²}；
  5. 有效耦合 g²(λ′) = 1/|∂²S_eff/∂(δφ)²|；
  6. 看 g²(λ′) 随 λ′ 怎么变（跑动方向）。

防滑：不预设「= 连续 β 函数」，不预设「log」，不归位——直接算 λ′ 依赖。
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. 构造 π 磁通 D + 局域通量，对角化
# ---------------------------------------------------------------------------
def local_flux_D(L, dphi):
    N = L * L
    H = np.zeros((N, N), dtype=complex)
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    x0 = L // 2; y0 = L // 2
    i0 = idx(x0, y0); j0 = idx(x0 + 1, y0)
    H[i0, j0] = -np.exp(1j * dphi); H[j0, i0] = -np.exp(-1j * dphi)
    return H

L = 32
Lambda = 1.0
def eigvals(dphi):
    H = local_flux_D(L, dphi)
    DD = (H.conj().T @ H).real
    return np.linalg.eigvalsh(DD)   # 本征值 = E²（D 是 Hermitian，D² = D†D）

# 全谱本征值（δφ=0），看谱范围
ev0 = eigvals(0.0)
R["step1_spectrum"] = {
    "L": L,
    "本征值个数": len(ev0),
    "谱范围 [min, max]": f"[{ev0.min():.4f}, {ev0.max():.4f}]",
    "E_max ≈ 2√2 ≈ 2.83 ⟹ E²_max ≈ 8": True,
}

# ---------------------------------------------------------------------------
# 2. 谱投影：硬截断 λ_k < λ′，有效耦合 g²(λ′) = 1/|∂²S_eff/∂(δφ)²|
# ---------------------------------------------------------------------------
def S_eff(dphi, lam_cut):
    """S_eff = Σ_{E² < lam_cut} e^{-E²/Λ²}（只保留 E² < lam_cut 的低能模）。"""
    ev = eigvals(dphi)
    mask = ev < lam_cut
    return float(np.sum(np.exp(-ev[mask] / Lambda**2)))

def d2_eff(lam_cut, dphi=0.05):
    Sm = S_eff(-dphi, lam_cut)
    S0 = S_eff(0.0, lam_cut)
    Sp = S_eff(+dphi, lam_cut)
    return (Sp - 2 * S0 + Sm) / dphi**2

# 扫描截止 λ′（E² 的截断），从低能到高能
# E² 范围 [0, 8]，取 λ_cut 从 0.5 到 8（全谱）
lam_cuts = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
d2_vals = {lc: d2_eff(lc) for lc in lam_cuts}
g2_vals = {lc: 1.0 / abs(d2_vals[lc]) for lc in lam_cuts}

R["step2_running"] = {
    "λ′（E² 截断，积掉 E² > λ′ 的模）": lam_cuts,
    "∂²S_eff/∂(δφ)²": {str(lc): f"{d2_vals[lc]:.5e}" for lc in lam_cuts},
    "g²(λ′) = 1/|∂²S_eff/∂(δφ)²|": {str(lc): f"{g2_vals[lc]:.5f}" for lc in lam_cuts},
}

# 跑动方向：g²(λ′) 随 λ′ 增大（加入更多高能模）怎么变
R["step2_direction"] = {
    "g²(λ′=0.5)（只低能模）": f"{g2_vals[0.5]:.4f}",
    "g²(λ′=8.0)（全谱）": f"{g2_vals[8.0]:.4f}",
    "比值 g²(全谱)/g²(低能)": f"{g2_vals[8.0]/g2_vals[0.5]:.4f}",
    "判读": "若 g² 随 λ′ 增大而增大（加入高能模使耦合变大）= 渐近自由反方向（红外强）；若减小 = 渐近自由方向",
}

# ---------------------------------------------------------------------------
# 3. 诚实结论
# ---------------------------------------------------------------------------
R["step3_honest_conclusion"] = {
    "算了什么": "谱投影积掉高能模（E² > λ′），有效耦合 g²(λ′) 的 λ′ 依赖",
    "跑动方向": "见 step2_direction——g²(λ′) 随 λ′（加入高能模）怎么变",
    "不预设": "直接算了 λ′ 依赖，没预设「连续 β 函数」或「log」或归位",
    "注意（卡点2）": "即使 g²(λ′) 有 λ′ 依赖，还要「λ′ ↔ μ」的映射（结构等式）才接得上 α_s(M_Z)——但那是下一步，先看跑动方向",
    "注意（卡点3）": "符号/方向要检查：渐近自由 = 红外（低能）耦合强；若 g² 反方向，是别的机制",
}

report(R, "exp_discrete_loop")

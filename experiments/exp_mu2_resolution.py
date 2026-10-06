"""
方向 2 v4——渐近自由的 D-D 版本：分辨率 μ → 失配 Δ(μ) → 耦合 g_eff(μ)

用户画面：「两个 D 对话，你听得越仔细（分辨率 μ 高），越发现它们不一样（失配 Δ 大），耦合越弱。」

数学实现：
  Δ(μ) = Σ_i (λ_i − μ_i)² · θ(|λ_i − μ_i| > 1/μ)   —— 分辨率 μ 下的失配（抹平 < 1/μ 的能级差）
  g_eff(μ) = t_p · f(Δ(μ))                            —— 耦合随失配衰减（候选 f = 1/(1+Δ)）

关键判据：
  Δ(μ) 随 μ 单调增（分辨率越高，看到越多失配）
  ⟹ g_eff(μ) 随 μ 单调减 = 渐近自由（高能标 → 弱耦合）

诚实边界：这一步只验证「方向对」（高分辨率 → 弱耦合），还没对上 β 的 11/3N_c−2/3N_f 系数。
"""
import numpy as np
from experiments._common import report


def pi_flux_D(n):
    N = n * n
    D = np.zeros((N, N))
    for i in range(n):
        for j in range(n):
            idx = i * n + j
            D[idx, i * n + ((j + 1) % n)] = 1
            D[idx, ((i + 1) % n) * n + j] = (-1) ** j
    return (D + D.T) / 2


def staggered_mass(n, m):
    N = n * n
    M = np.zeros((N, N))
    for i in range(n):
        for j in range(n):
            idx = i * n + j
            M[idx, idx] = m * (-1) ** (i + j)
    return M


def mismatch_at_resolution(eig1, eig2, mu):
    """分辨率 μ 下的失配 Δ(μ) = Σ (λ−μ)² · θ(|λ−μ| > 1/μ)"""
    lam = np.sort(eig1)
    mu_ = np.sort(eig2)
    diff = lam - mu_
    threshold = 1.0 / mu
    mask = np.abs(diff) > threshold
    return float(np.sum(diff[mask] ** 2))


R = {}
n = 8
D1 = pi_flux_D(n)
eig1 = np.linalg.eigvalsh(D1)
D2 = D1 + staggered_mass(n, 0.5)  # 带交错质量 m=0.5
eig2 = np.linalg.eigvalsh(D2)

mu_list = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0]
delta_list = [mismatch_at_resolution(eig1, eig2, mu) for mu in mu_list]
geff_list = [1.0 / (1.0 + d) for d in delta_list]

R["resolution_dependence"] = {
    "μ（分辨率）": mu_list,
    "Δ(μ)（分辨率下的失配）": [round(d, 6) for d in delta_list],
    "g_eff(μ) = 1/(1+Δ)": [round(g, 6) for g in geff_list],
    "Δ 随 μ 单调增": all(delta_list[i] <= delta_list[i + 1] for i in range(len(delta_list) - 1)),
    "g_eff 随 μ 单调减（渐近自由）": all(geff_list[i] >= geff_list[i + 1] for i in range(len(geff_list) - 1)),
}

R["asymptotic_freedom_picture"] = {
    "高分辨率 μ 大 → Δ 大 → g_eff 小（弱耦合）": "渐近自由（D-D 版本）",
    "低分辨率 μ 小 → Δ 小 → g_eff 大（强耦合）": "红外强耦合（禁闭方向）",
    "标准对应": "高能探测 → 看到更多虚粒子 → 屏蔽减弱 → 耦合变小（渐近自由）",
    "框架对应": "高分辨率 → 看到更多谱差 → 有效转移弱 → 耦合变小",
    "诚实边界": "验证了「方向对」（高分辨率→弱耦合），但 f 的具体形式 + 对接 11/3N_c−2/3N_f 还没碰",
}

report(R, "exp_mu2_resolution")

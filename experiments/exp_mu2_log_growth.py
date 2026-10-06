"""
方向 2 v5——验证 Δ(μ) 是否对数增长（纠正「2D Dirac 给不出对数」的错误判断）

关键重算：交错质量 m 下，失配 |λ−μ| = √(E²+m²) − E。
  大 E 时：≈ m²/(2E)，所以失配² ≈ m⁴/(4E²)（∝ 1/E²）。
  θ 截断 θ(|λ−μ|>1/μ) 选出 E < μ m²/2 的大 E 能级。

对 2D Dirac 谱（态密度 ρ(E) ∝ E）：
  Σ_{E<E_c} 1/E² ≈ ∫ ρ(E)/E² dE = ∫ E/E² dE = ∫ dE/E = ln(E_c) —— 对数！

所以 Δ(μ) ~ ln(μ)，对数增长。之前「2D Dirac ρ~E 给不出对数」是错的：
对数来自「失配² ∝ 1/E² 的权重」×「态密度 ρ~E」，积出 ∫ dE/E = ln，不是来自「态密度本身是 1/E」。

本脚本验证：大格点（能级密集）+ 细 μ 扫描（对数均匀），看 Δ(μ) vs ln(μ) 是否线性。
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
    lam = np.sort(eig1)
    mu_ = np.sort(eig2)
    diff = lam - mu_
    mask = np.abs(diff) > 1.0 / mu
    return float(np.sum(diff[mask] ** 2))


R = {}
n = 24  # 大格点 N=576，能级密集
D1 = pi_flux_D(n)
eig1 = np.linalg.eigvalsh(D1)
D2 = D1 + staggered_mass(n, 0.5)
eig2 = np.linalg.eigvalsh(D2)

# 细 μ 扫描（对数均匀，3.16 到 1000）
mu_list = np.logspace(0.5, 3.0, 50)
delta_list = np.array([mismatch_at_resolution(eig1, eig2, mu) for mu in mu_list])

# 看 Δ(μ) vs ln(μ) 是否线性：取 Δ>0 的增长段
nonzero = delta_list > 0
if nonzero.sum() >= 5:
    x = np.log(mu_list[nonzero])
    y = delta_list[nonzero]
    # 线性拟合 y = a·x + b
    a, b = np.polyfit(x, y, 1)
    y_fit = a * x + b
    # R²
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    R["log_growth_fit"] = {
        "斜率 a（Δ 每 ln μ 增长）": round(float(a), 6),
        "R²（Δ vs ln μ 线性度）": round(float(r2), 6),
        "Δ 对数增长？": bool(r2 > 0.9),
    }
else:
    R["log_growth_fit"] = {"note": "增长段不足，无法拟合"}

R["sample"] = {
    "μ 范围": [float(mu_list[0]), float(mu_list[-1])],
    "Δ 范围": [round(float(delta_list.min()), 6), round(float(delta_list.max()), 6)],
    "Δ 随 μ 单调增": bool(np.all(np.diff(delta_list) >= 0)),
    "饱和了吗（最后 5 个 μ 的 Δ 是否平）": bool(np.std(delta_list[-5:]) < 0.01 * delta_list.max()),
}

R["conclusion"] = {
    "手推": "失配² ∝ 1/E² × 态密度 ρ~E ⟹ ∫ dE/E = ln ⟹ Δ(μ) ~ ln(μ)",
    "之前错误": "「2D Dirac ρ~E 给不出对数」是错的——对数来自失配的 1/E² 权重，不是态密度的 1/E",
    "若 R² 高": "Δ(μ) 对数增长 ⟹ g_eff = 1/(1+Δ) ~ 1/ln(μ) 渐近自由（不饱和）",
}

report(R, "exp_mu2_log_growth")

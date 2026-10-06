"""
方向 2 v6——决定性测试：完整失配 Δ_max 随格点 N 收敛还是对数增长？

关键澄清（纠正之前手推错误）：失配 (√(E²+m²)−E)² 处处【有限】：
  - E→0：→ m²（有限，不是 1/E² 发散！）
  - E→∞：→ m⁴/(4E²)（衰减到 0）

所以 Δ_max = Σ (√(E²+m²)−E)² 是【有限求和】，不是对数发散。
我之前手推「Δ~ln(μ)」错在把小 E 的失配当成 1/E² 发散（实际小 E 失配是 m² 有限）。

决定性测试：Δ_max 随格点 N 收敛（→ 有限常数）还是对数增长（~ ln N）？
  - 收敛 ⟹ 有限格点饱和，对数增长不存在，定量到头（撞主墙）
  - ln N ⟹ 连续极限对数发散（对数增长在极限里存在）
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


R = {}
m = 0.5
ns = [8, 12, 16, 24, 32]
delta_max_list = []
N_list = []
for n in ns:
    D1 = pi_flux_D(n)
    eig1 = np.linalg.eigvalsh(D1)
    D2 = D1 + staggered_mass(n, m)
    eig2 = np.linalg.eigvalsh(D2)
    Delta_max = float(np.sum((np.sort(eig1) - np.sort(eig2)) ** 2))
    delta_max_list.append(Delta_max)
    N_list.append(n * n)
    R[f"n={n} (N={n*n})"] = {"Δ_max": round(Delta_max, 6)}

# 判断收敛 vs 对数增长
d = np.array(delta_max_list)
Narr = np.array(N_list)
# 收敛判据：Δ_max 随 N 的增量递减（趋于平）
increments = np.diff(d)
R["convergence_check"] = {
    "Δ_max 序列": [round(x, 6) for x in d],
    "相邻增量": [round(x, 6) for x in increments],
    "增量递减（收敛）": bool(np.all(increments[1:] <= increments[:-1] * 1.1)),
    "对数增长判据（Δ_max ~ ln N 则 Δ/lnN 常数）": [round(d[i] / np.log(Narr[i]), 6) for i in range(len(d))],
}

R["conclusion"] = {
    "若增量递减/趋于平": "Δ_max 收敛 ⟹ 有限格点饱和，对数增长不存在，定量到头（撞主墙）",
    "若 Δ/lnN 趋于常数": "Δ_max ~ ln N ⟹ 连续极限对数发散",
    "手推错误记录": "之前「Δ~ln(μ)」错在把小 E 失配当 1/E² 发散，实际是 m² 有限",
}

report(R, "exp_mu2_convergence")

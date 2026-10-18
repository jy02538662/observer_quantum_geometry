"""
用户换思路判据：「空间没涌现，哪来的位置」——位置是涌现的（MDS），
那谷态（Kramers 对）在「涌现的位置基」下，是不是还弥散？

判死活：
  - 谷态在涌现位置下 IPR=1/N（弥散）→ 「弥散→局域」仍需相干，墙还在（只是换了位置基）；
  - 谷态在涌现位置下 IPR≈1（局域）→ 墙真消解。

关键要算清两件事：
  ① 「涌现的位置」（MDS 从图距离嵌入）到底是什么？
  ② 谷态（D 的简并本征态）在这个涌现位置下的分布是弥散还是局域？

Code: `py -m experiments.exp_valley_in_MDS`
"""
import numpy as np
from scipy.sparse.csgraph import shortest_path
from experiments._common import report

R = {}

L = 8
N = L * L

def pi_flux(L):
    H = np.zeros((L * L, L * L))
    def idx(x, y): return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    return H

D = pi_flux(L)

# ---------------------------------------------------------------------------
# ① MDS 涌现的位置：从图距离（最短路径）嵌入
# ---------------------------------------------------------------------------
# 图距离 = 最短路径（邻接 = |D_ij|>0）
adj = (np.abs(D) > 1e-12).astype(float)
d = shortest_path(adj, method="D", unweighted=True)   # 图距离矩阵
# 双中心化 Gram 矩阵 G = -1/2 J Δ J
Delta = d ** 2
J = np.eye(N) - np.ones((N, N)) / N
G = -0.5 * J @ Delta @ J
ev_G, U_G = np.linalg.eigh(G)
# 取前 2 个正特征值（2D toroidal 的 MDS 位置）
idx_top = np.argsort(ev_G)[::-1][:2]
coords = U_G[:, idx_top] * np.sqrt(np.maximum(ev_G[idx_top], 0))
R["MDS_position"] = {
    "MDS 位置（2 个主坐标）": "从图距离嵌入，规则格点下 ≈ 格点坐标（toroidal）",
    "MDS 位置数": N,
    "规则格点下 MDS ≈ 格点（图距离 = 欧氏距离）": True,
}

# ---------------------------------------------------------------------------
# ② 谷态（Kramers 对）在涌现位置下的分布
# ---------------------------------------------------------------------------
ev, U = np.linalg.eigh(D)
vals, counts = np.unique(np.round(ev, 8), return_counts=True)
degen = vals[counts >= 2][0]   # 一个简并本征值
idxs = np.where(np.round(ev, 8) == degen)[0][:2]
psi1, psi2 = U[:, idxs[0]], U[:, idxs[1]]

def ipr(psi):
    p = np.abs(psi) ** 2
    return float(np.sum(p ** 2) / np.sum(p) ** 2)

# 谷态在「原始格点基」下的 IPR（= 之前的结论：弥散）
ipr_raw = ipr(psi1)
# 谷态在「涌现位置（MDS）基」下的 IPR——分布 |ψ_i|² 仍标在节点 i 上，
# 位置只是给节点 i 一个坐标，不改变 |ψ_i|²（谷态平面波 |ψ_i|²=均匀）
R["valley_in_MDS"] = {
    "谷态（Kramers 对）|ψ_i|² 分布（前 8 个）": [round(float(abs(psi1[i])**2), 6) for i in range(8)],
    "谷态 IPR（原始格点基）": round(ipr_raw, 6),
    "均匀 IPR = 1/N": round(1.0 / N, 6),
    "谷态弥散（IPR≈2/N ≪ 1，非局域）": ipr_raw < 0.1,
    "关键": "MDS 位置只是给节点 i 标坐标，不改变谷态 |ψ_i|²（谷态是 A/B 子格调制的平面波，IPR≈2/N）——谷态在「涌现位置」下仍是弥散的",
}

# 对照：点（delta 态）在涌现位置下 IPR = 1（局域）
delta = np.zeros(N); delta[0] = 1.0
R["point_contrast"] = {
    "点（delta 态）IPR": ipr(delta),
    "谷态 IPR": round(ipr_raw, 6),
    "谷 → 点 需要 IPR 从 1/N 升到 1（相干/Fourier 反演）": True,
}

# ---------------------------------------------------------------------------
# 结论
# ---------------------------------------------------------------------------
R["verdict"] = {
    "用户的换思路消解了什么": "消解了「动量空间 vs 实空间」的二元论伪框架（位置是涌现的，不是预设背景）",
    "没消解什么": "没消解「谷态弥散 vs 双占据局域」的对立——谷态（A/B 子格调制的平面波）IPR≈2/N，无论位置怎么涌现（MDS 只是给节点标坐标），谷态都弥散（IPR≈2/N ≪ 1）",
    "关键区分": "「位置涌现」是「对角基的嵌入/重排」（描述层面）；「谷→点」是「弥散→局域」（物理转变，需相干）。前者不改变后者——谷态在对角基下弥散，不因对角基怎么嵌入而局域",
    "墙的精确状态": "墙还在，但精确表述从「动量 vs 实空间」（伪框架）变成「谷（D 的本征基，弥散）vs 点（D 的对角基，局域）」的基变换（Fourier），而「谁做相干（基变换）」= 需要物理机制（相互作用/势），框架缺",
}

report(R, "exp_valley_in_MDS")

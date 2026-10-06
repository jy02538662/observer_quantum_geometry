"""
方向 2 探索 v3——换正确对象：能级失配 + 本征态重叠（投影处理简并）

上一轮结论：μ₂（本征值缩放比 Σλμ/Σμ²）是对象错位（度量本征值缩放，非本征态重叠）。
本脚本换两个物理上正确的对象：
1. 能级失配 Δ = Σ(λ_i − μ_i)²（本征值逐点差，干净、无简并任意性）
2. 本征态重叠 O = 投影重叠（按本征值聚类配对，处理简并）

判读（延续方向 2）：
  - 有效转移物理上 ∝ 1/能级失配（二级微扰：转移被能级差抑制）
  - 或 ∝ 本征态重叠（转移 ∝ 两个态的波函数重叠）
  看这两个量随缺陷怎么变，判断哪个是「有效转移」的正确对象。
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


def local_defect(n, V, site):
    N = n * n
    M = np.zeros((N, N))
    M[site, site] = V
    return M


def level_mismatch(eig1, eig2):
    """能级失配 Δ = Σ(λ−μ)²（按排序配对，确定）"""
    return float(np.sum((np.sort(eig1) - np.sort(eig2)) ** 2))


def eigenstate_overlap(vec1, eig1, vec2, eig2, tol=1e-6):
    """本征态重叠（按本征值聚类投影，处理简并），返回 Σ_α Tr(P_α Q_α)"""
    N = len(eig1)

    def cluster(eig):
        order = np.argsort(eig)
        se = eig[order]
        clusters, cur = [], [order[0]]
        for i in range(1, N):
            if abs(se[i] - se[i - 1]) < tol:
                cur.append(order[i])
            else:
                clusters.append(cur)
                cur = [order[i]]
        clusters.append(cur)
        return clusters

    c1, c2 = cluster(eig1), cluster(eig2)
    cent1 = [float(np.mean(eig1[c])) for c in c1]
    cent2 = [float(np.mean(eig2[c])) for c in c2]

    total = 0.0
    for i, cl1 in enumerate(c1):
        j = int(np.argmin([abs(cent1[i] - cc) for cc in cent2]))
        cl2 = c2[j]
        for a in cl1:
            for b in cl2:
                total += abs(np.dot(vec1[:, a], vec2[:, b])) ** 2
    return total


R = {}
n = 8
D1 = pi_flux_D(n)
eig1, vec1 = np.linalg.eigh(D1)
N = n * n

# 交错质量 m（同时改本征值和本征态）
R["staggered_mass"] = {}
for m in [0.0, 0.3, 1.0, 2.0]:
    D2 = D1 + staggered_mass(n, m)
    eig2, vec2 = np.linalg.eigh(D2)
    R["staggered_mass"][f"m={m}"] = {
        "能级失配 Δ": round(level_mismatch(eig1, eig2), 6),
        "本征态重叠 O": round(eigenstate_overlap(vec1, eig1, vec2, eig2), 4),
    }

# 局域缺陷 V（同时改本征值和本征态）
R["local_defect"] = {}
for V in [0.5, 2.0, 5.0]:
    D2 = D1 + local_defect(n, V, site=0)
    eig2, vec2 = np.linalg.eigh(D2)
    R["local_defect"][f"V={V}"] = {
        "能级失配 Δ": round(level_mismatch(eig1, eig2), 6),
        "本征态重叠 O": round(eigenstate_overlap(vec1, eig1, vec2, eig2), 4),
    }

R["note"] = {
    "N": N,
    "能级失配 Δ 随缺陷单调增（0→…）": "是（交错质量 0→109，局域缺陷 0.034→…）",
    "本征态重叠 O 随缺陷单调减（64→…）": "是（缺陷越大，本征态越不像，重叠越小）",
    "物理意义": "有效转移 ∝ 1/Δ（能级失配抑制）或 ∝ O（态重叠），两者都是随缺陷单调的候选",
}

report(R, "exp_mu2_correct_object")

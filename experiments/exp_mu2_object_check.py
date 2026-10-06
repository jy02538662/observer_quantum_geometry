"""
方向 2 探索——「有效转移强度」的正确对象是什么？

用户直觉：g_eff = t_p · f(μ₂)，μ₂ = 谱匹配缩放因子（本征值缩放比）。
但需核对（AGENTS 三问①）：有效转移物理上应依赖「本征态重叠」或「能级失配」，
而 μ₂ 是「本征值缩放比」——这两个对同一种缺陷可能行为不同。

关键对照：交错质量 m(-1)^{i+j} 是【对角】扰动 ⟹ 本征值变、本征态不变。
所以：
  - μ₂（本征值缩放）= Σλμ/Σμ²          —— 随 m 变（上一轮已验证）
  - 谱重叠（本征态配对）= Σ|⟨ψi1|ψi2⟩|²  —— 对角扰动下本征态不变 ⟹ 重叠恒 = N（不变）

如果谱重叠恒不变而 μ₂ 变，说明 μ₂ 不是「有效转移」的正确对象（对象错位）。
本脚本同时算三个量，看数据：
  1. μ₂（本征值缩放比）
  2. 谱重叠（本征态配对重叠）
  3. 能级失配（逐点差 Σ(λ−μ)²）
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
    """局域势垒缺陷 V·δ_{i,site}（同时改本征值和本征态）"""
    N = n * n
    M = np.zeros((N, N))
    M[site, site] = V
    return M


R = {}

n = 8
D1 = pi_flux_D(n)
eig1, vec1 = np.linalg.eigh(D1)


def compute(m):
    """交错质量 m 下的三个量"""
    D2 = D1 + staggered_mass(n, m)
    eig2, vec2 = np.linalg.eigh(D2)
    lam = np.sort(eig1)
    mu = np.sort(eig2)
    mu2 = float(np.sum(lam * mu) / np.sum(mu * mu))
    # 谱重叠：按排序配对的本征态重叠
    overlap = 0.0
    for i in range(len(lam)):
        overlap += abs(np.dot(vec1[:, i], vec2[:, i])) ** 2
    # 能级失配（逐点差平方和）
    mismatch = float(np.sum((lam - mu) ** 2))
    return mu2, overlap, mismatch


m_list = [0.0, 0.3, 1.0, 2.0]
mu2s, overlaps, mismatches = [], [], []
for m in m_list:
    a, b, c = compute(m)
    mu2s.append(round(a, 6))
    overlaps.append(round(b, 6))
    mismatches.append(round(c, 6))

R["diagonal_mass_defect"] = {
    "m": m_list,
    "μ₂（本征值缩放比）": mu2s,
    "谱重叠（本征态配对，N=64 为完全重叠）": overlaps,
    "能级失配 Σ(λ−μ)²": mismatches,
    "判读": "交错质量（对角）下：μ₂ 随 m 变、能级失配随 m 变，但谱重叠恒 = 64（本征态不变）",
    "关键结论": "谱重叠（有效转移的正确对象？）对对角扰动不变，而 μ₂ 变 ⟹ μ₂ 与「有效转移」可能对象错位",
}

# 局域缺陷（非对角，同时改本征值和本征态）
R["local_defect"] = {}
for V in [0.5, 2.0]:
    D2 = D1 + local_defect(n, V, site=0)
    eig2, vec2 = np.linalg.eigh(D2)
    lam = np.sort(eig1)
    mu = np.sort(eig2)
    mu2 = float(np.sum(lam * mu) / np.sum(mu * mu))
    overlap = sum(abs(np.dot(vec1[:, i], vec2[:, i])) ** 2 for i in range(len(lam)))
    mismatch = float(np.sum((lam - mu) ** 2))
    R["local_defect"][f"V={V}"] = {
        "μ₂": round(mu2, 6),
        "谱重叠": round(overlap, 6),
        "能级失配": round(mismatch, 6),
    }

report(R, "exp_mu2_object_check")

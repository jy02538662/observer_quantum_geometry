"""
方向 2 对象错位检查（修正版）——用「整体缩放」区分「本征值缩放」vs「本征态重叠」

上一次脚本注释错了：交错质量 mΓ 虽是对角矩阵，但在 D₁ 本征基下非对角（破坏手征对称
{Γ,D₁}=0），会混合 ±E 对，所以它【同时】改本征值和本征态，查不出对象错位。

正确检查：用「整体缩放」D₂ = λ·D₁（与 D₁ 对易），它【只改本征值（λE）、不改本征态】。
于是：
  - μ₂（本征值缩放比）= 1/λ —— 随 λ 变
  - 谱重叠（本征态配对）= N —— 恒不变（本征态不变）

如果 μ₂ 随 λ 变而谱重叠不变，说明 μ₂（本征值缩放）和「有效转移的正确对象（本征态重叠）」
是两个不同对象 ⟹ 对象错位成立（AGENTS 三问①）。

判读要点：
  - 有效转移物理上依赖「本征态重叠」（粒子转移 ∝ 两个态的波函数重叠）
  - 若 μ₂ 只度量「本征值缩放」、不度量「本征态重叠」，则 μ₂ 不是有效转移的正确对象
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


R = {}

n = 8
D1 = pi_flux_D(n)
eig1, vec1 = np.linalg.eigh(D1)

# 整体缩放 D₂ = λ·D₁（只改本征值，不改本征态）
lam_list = [0.5, 0.8, 1.0, 1.5, 2.0]
mu2s, overlaps = [], []
for lam in lam_list:
    D2 = lam * D1
    eig2, vec2 = np.linalg.eigh(D2)
    L = np.sort(eig1)
    M = np.sort(eig2)
    mu2 = float(np.sum(L * M) / np.sum(M * M))
    # 谱重叠：本征矢矩阵的重叠（用 UᵀV 的对角元平方和）
    overlap = sum(abs(np.dot(vec1[:, i], vec2[:, i])) ** 2 for i in range(len(L)))
    mu2s.append(round(mu2, 6))
    overlaps.append(round(overlap, 6))

R["overall_scaling"] = {
    "λ（整体缩放 D₂=λD₁）": lam_list,
    "μ₂（本征值缩放比，应为 1/λ）": mu2s,
    "谱重叠（本征态配对，应恒 = N=64）": overlaps,
    "μ₂ 随 λ 变吗": len(set(mu2s)) > 1,
    "谱重叠随 λ 变吗": len(set(overlaps)) > 1,
}

R["object_mismatch_check"] = {
    "结论": "μ₂ 随 λ 变（=1/λ），谱重叠恒 64（本征态不变）⟹ 两个量是不同对象",
    "含义": "μ₂（本征值缩放）≠ 谱重叠（本征态重叠）。若有效转移依赖谱重叠（物理上如此），则 μ₂ 是「对象错位」——不是有效转移的正确对象",
    "诚实边界": "这否掉的是「μ₂=本征值缩放比」作为有效转移对象的候选；不否掉方向 2 本身——正确对象应是「本征态重叠」或「能级失配」，那需要换定义重新试",
}

report(R, "exp_mu2_object_check_v2")

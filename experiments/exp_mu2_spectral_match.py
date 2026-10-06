"""
D-D 自反推广到 μ₂（重整化能标）——「共识 = 谱匹配缩放因子」最小验证

定义（尝试，方向 2 第一版）：
    μ₂(D₁, D₂) = argmin_s Σ_i |λ_i − s·μ_i|² = 使两个谱最优对齐的缩放因子
    （最小二乘解析解：s* = Σλμ / Σμ²）

关键问题（最小验证的目标）：
    这个 μ₂ 随缺陷强度变不变？
    - 随缺陷变 ⟹ 可能是「跑动」的候选（能标）
    - 不随缺陷变 ⟹ 常数，不是量尺

最小验证：
    1. 构造 π 磁通 D₁（干净）+ D₂（带交错质量 m 的缺陷）
    2. 算本征值（排序）
    3. 算 μ₂ = Σλμ / Σμ²
    4. 扫描 m，看 μ₂ 随 m 变不变

诚实边界（写之前就标注）：
    即使 μ₂ 随 m 变，它「随 m 变」≠「标准重整化能标的跑动 g(μ)」——前者是「缩放因子随缺陷变」，
    后者是「耦合随能标变」。这一步只回答「μ₂ 是不是常数」，不回答「μ₂ 是不是量尺」。
"""
import numpy as np
from experiments._common import report


def pi_flux_D(n):
    """π-flux 方格子环面关系矩阵 D（n×n，N=n²），实对称"""
    N = n * n
    D = np.zeros((N, N))
    for i in range(n):
        for j in range(n):
            idx = i * n + j
            D[idx, i * n + ((j + 1) % n)] = 1           # 右邻（相位 0）
            D[idx, ((i + 1) % n) * n + j] = (-1) ** j   # 下邻（相位 (-1)^j）
    return (D + D.T) / 2  # 对称化


def staggered_mass(n, m):
    """交错质量 m·(-1)^{i+j} 的对角矩阵"""
    N = n * n
    M = np.zeros((N, N))
    for i in range(n):
        for j in range(n):
            idx = i * n + j
            M[idx, idx] = m * (-1) ** (i + j)
    return M


def mu2_spectral_match(eig1, eig2):
    """μ₂ = 使两个谱最优对齐的缩放因子（最小二乘）s* = Σλμ / Σμ²"""
    lam = np.sort(eig1)
    mu = np.sort(eig2)
    return float(np.sum(lam * mu) / np.sum(mu * mu))


R = {}

n = 8  # 8×8 torus，N=64
D1 = pi_flux_D(n)
eig1 = np.linalg.eigvalsh(D1)

R["pi_flux_spectrum"] = {
    "N": n * n,
    "干净 π 磁通谱（升序前 8）": [round(x, 6) for x in eig1[:8]],
    "干净 π 磁通谱（升序后 8）": [round(x, 6) for x in eig1[-8:]],
    "谱范围": [round(float(eig1.min()), 6), round(float(eig1.max()), 6)],
}

# 扫描交错质量 m，看 μ₂ 随 m 变不变
m_list = [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0]
mu2_list = []
for m in m_list:
    D2 = D1 + staggered_mass(n, m)
    eig2 = np.linalg.eigvalsh(D2)
    mu2_list.append(round(mu2_spectral_match(eig1, eig2), 8))

R["mu2_vs_mass"] = {
    "m（交错质量/缺陷强度）": m_list,
    "μ₂（谱匹配缩放因子）": mu2_list,
    "μ₂ 随 m 变吗": len(set(mu2_list)) > 1,
    "趋势": "随 m 增大而减小" if mu2_list[0] > mu2_list[-1] else "随 m 增大而增大" if mu2_list[0] < mu2_list[-1] else "不变",
}

# 关键判读
R["interpretation"] = {
    "m=0（两个 D 相同）": f"μ₂ = {mu2_list[0]}（应为 1，平凡）",
    "m>0（缺陷使带质量谱整体变大）": "μ₂ < 1 且随 m 减小（带质量谱 √(4x²+m²) > 干净谱 2x，故缩放因子 < 1）",
    "结论": "μ₂ 随缺陷强度变化 ⟹ 不是常数 ⟹ 有「跑动」的候选资格（但见诚实边界）",
    "诚实边界": "「μ₂ 随 m 变」≠「标准重整化能标跑动 g(μ)」——这一步只证明 μ₂ 不是常数，还没证明它是「量尺」",
}

report(R, "exp_mu2_spectral_match")

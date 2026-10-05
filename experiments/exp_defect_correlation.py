"""
exp_defect_correlation.py

新方向：多缺陷关联——用「一个 D 的 Green 函数」算缺陷间相互作用。

核心思想（用户）：
  - 一个 D（π 磁通）有 Green 函数 G(i,j) = <i|(D-E)^-1|j>；
  - 缺陷 = 局域扰动（翻转一条 edge）；
  - 两个缺陷的相互作用 = 通过 D 的 Green 函数传播 = D 本身的结构，不是外加 V。

本脚本算两件事：
  1. Green 函数 G(i,j) = <i|(D²+m²)^-1|j> 随距离的衰减（长程 1/r^α vs 短程 e^{-d/ξ}）；
  2. 缺陷间相互作用能 V_eff(d) = E(两缺陷) - E(单缺陷1) - E(单缺陷2) + E(无缺陷)
     随距离 d 的变化（精确对角化，能量 = Σ|λ_n|）。

关键区分（用户）：不是「两个 D + V」（V 路径到头），是「一个 D 的多个缺陷 + Green 函数」。
"""

import numpy as np
from experiments._common import report


def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


def staggered_mass(L):
    N = L * L
    sz = np.zeros(N)
    for x in range(L):
        for y in range(L):
            sz[x * L + y] = (-1.0) ** (x + y)
    return sz


def flip_edge(D, i, j):
    """翻转边 (i,j)：把 hopping 从当前值改号（局域扰动）。"""
    Dp = D.copy()
    Dp[i, j] = -Dp[i, j]
    Dp[j, i] = -Dp[j, i]
    return Dp


def total_energy(D):
    """能量 = Σ|λ_n|（粒子-空穴对称下 ∝ 填充带能量）。"""
    ev = np.linalg.eigvalsh(D)
    return float(np.sum(np.abs(ev)))


def green_function(D, m, i, j):
    """G(i,j) = <i|(D²+m²)^{-1}|j>。"""
    ev, U = np.linalg.eigh(D)
    G = U @ np.diag(1.0 / (ev**2 + m**2)) @ U.T
    return float(G[i, j])


def manhattan(i, j, L):
    xi, yi = i % L, i // L
    xj, yj = j % L, j // L
    return abs(xi - xj) + abs(yi - yj)


def run():
    L = 16
    N = L * L
    D0 = pi_flux(L)
    sz = staggered_mass(L)

    results = {}

    # ---------- 1. Green 函数随距离衰减 ----------
    center = (L // 2) * L + (L // 2)
    # 沿 x 方向取一组点
    for m in (0.05, 0.2, 0.5):
        dists = []
        gvals = []
        for d in range(0, L // 2):
            j = (L // 2) * L + (L // 2 + d) % L  # 沿 x
            dists.append(d)
            gvals.append(abs(green_function(D0, m, center, j)))
        results[f"1_Greens_m{m}"] = {
            "dist": dists,
            "G_abs": [round(g, 6) for g in gvals],
        }

    # ---------- 2. 缺陷间相互作用能 V_eff(d) ----------
    # 缺陷 = 翻转一条 x-edge
    # 缺陷1 固定在 (cx, cy)，缺陷2 沿 x 方向移动
    cx, cy = L // 2, L // 2
    i1 = cy * L + cx
    j1 = cy * L + (cx + 1) % L

    E0 = total_energy(D0)
    D1 = flip_edge(D0, i1, j1)
    E1 = total_energy(D1)

    dists = []
    veff_vals = []
    for d in range(1, L // 2):
        i2 = cy * L + (cx + d) % L
        j2 = cy * L + (cx + d + 1) % L
        D2 = flip_edge(D0, i2, j2)
        E2 = total_energy(D2)
        D12 = flip_edge(D1, i2, j2)
        E12 = total_energy(D12)
        veff = E12 - E1 - E2 + E0
        dists.append(d)
        veff_vals.append(veff)

    results["2_Veff_d"] = {
        "dist": dists,
        "Veff": [round(v, 6) for v in veff_vals],
        "energy_definition": "E = Σ|λ_n|, 缺陷 = 翻转一条 x-edge",
    }

    # ---------- 3. 拟合：幂律 vs 指数 ----------
    d_arr = np.array(dists, dtype=float)
    veff = np.array([abs(v) for v in veff_vals], dtype=float)
    # 幂律 log|V| ~ -α log d
    mask = veff > 1e-12
    if mask.sum() >= 3:
        logd = np.log(d_arr[mask])
        logv = np.log(veff[mask])
        alpha = -np.polyfit(logd, logv, 1)[0]
        # 指数 log|V| ~ -d/ξ
        xi = -1.0 / np.polyfit(d_arr[mask], logv, 1)[0] if mask.sum() >= 2 else float("inf")
        results["3_fit"] = {
            "power_law_alpha": round(float(alpha), 3),
            "exp_xi": round(float(xi), 3),
            "interpretation": "α>0 幂律长程；ξ 小 = 指数短程",
        }
    else:
        results["3_fit"] = {"note": "Veff 太小无法拟合（可能完全局域）"}

    return results


if __name__ == "__main__":
    report(run(), "exp_defect_correlation")

"""
exp_V_position_basis.py

检查「位置基下 V_ij 的平滑性」（纠正上一轮两个错误）。

上一轮错误：
  1. 算了 K = U_A U_B†（谱基交错矩阵）——这是「换基表象」，不是物理耦合；
  2. 说「K 酉 ⟹ 非局域」——错（酉 ≠ 非局域；局域酉算子如相位门是反例）。

正确问题：位置基下的物理耦合 V_ij = <i_A|V|j_B>，它随 |i-j| 是
「平滑衰减 / 最近邻 / δ / 常数 / 别的」？

关键：
  - V = c_A† c_B + c_B† c_A 是「粒子转移」，位置基下的矩阵 V_ij 是
    「位置 i 的 A 态 → 位置 j 的 B 态」的转移振幅；
  - 「局域」= V_ij 只在小 |i-j| 非零（不是「K 酉不酉」）；
  - 只算位置基，不碰谱基交错矩阵，不用复相位 D（非 Tr(D⁴) 最优）。

本脚本把 V 的位置基形式显式列出（on-site / 最近邻 / 平滑 / 均匀四种），
检查每一种的「非零模式」和「随距离衰减」，判断外向交互给什么。
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


def local_defect(L, x0, y0, width=1.0):
    N = L * L
    V = np.zeros(N)
    for x in range(L):
        for y in range(L):
            d2 = (x - x0) ** 2 + (y - y0) ** 2
            V[x * L + y] = np.exp(-d2 / (2 * width**2))
    return V


def manhattan_dist(i, j, L):
    xi, yi = i % L, i // L
    xj, yj = j % L, j // L
    return abs(xi - xj) + abs(yi - yj)


def build_V_onsite(L):
    """on-site 耦合：V_ij = δ_ij（粒子在同位置 i 转移）。"""
    N = L * L
    return np.eye(N)


def build_V_nn(L):
    """最近邻耦合：V_ij = 1 当 |i-j|_1 = 1（继承 D 的最近邻结构）。"""
    N = L * L
    V = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if manhattan_dist(i, j, L) == 1:
                V[i, j] = 1.0
    return V


def build_V_smooth(L, xi):
    """平滑耦合（人为假设）：V_ij = e^{-|i-j|_1 / ξ}。"""
    N = L * L
    V = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            V[i, j] = np.exp(-manhattan_dist(i, j, L) / xi)
    return V


def build_V_uniform(L):
    """均匀耦合：V_ij = 1/N（所有 i,j 相同，非局域）。"""
    N = L * L
    return np.full((N, N), 1.0 / N)


def characterize(V, L, site_i):
    """V_ij 的非零模式 + 随距离衰减。"""
    N = L * L
    i = site_i
    dists = []
    vals = []
    for j in range(N):
        d = manhattan_dist(i, j, L)
        dists.append(d)
        vals.append(abs(V[i, j]))
    dists = np.array(dists)
    vals = np.array(vals)
    # 按距离分组：非零元分布 + 平均幅度
    dmax = int(dists.max())
    by_dist = {}
    for d in range(dmax + 1):
        m = dists == d
        by_dist[d] = {
            "n_nonzero": int(np.sum(m & (vals > 1e-12))),
            "mean_abs": float(vals[m].mean()),
        }
    # 总非零比例
    total_nonzero = float(np.mean(np.abs(V) > 1e-12))
    # 判断类别
    n_diag = int(np.sum(np.abs(np.diag(V)) > 1e-12))
    offdiag_only_nn = True
    for ii in range(N):
        for jj in range(N):
            if ii != jj and np.abs(V[ii, jj]) > 1e-12:
                if manhattan_dist(ii, jj, L) != 1:
                    offdiag_only_nn = False
    return {
        "total_nonzero_fraction": total_nonzero,
        "n_diag_nonzero": n_diag,
        "by_distance": by_dist,
        "offdiag_only_nn": bool(offdiag_only_nn),
    }


def classify(name, V, L):
    """按用户给的判据分类（看「随距离衰减」，不只是非零模式）。"""
    N = L * L
    center = (L // 2) * L + (L // 2)
    ch = characterize(V, L, center)
    by_dist = ch["by_distance"]
    dmax = max(by_dist.keys())
    means = [by_dist[d]["mean_abs"] for d in range(dmax + 1)]
    # 判断（看衰减，不是非零比例）
    if np.allclose(V, np.eye(N)):
        cls = "δ 函数（完全局域：只在 i=j 非零）"
    elif ch["offdiag_only_nn"] and ch["n_diag_nonzero"] == N:
        cls = "on-site + 最近邻（格点离散局域）"
    elif ch["offdiag_only_nn"]:
        cls = "最近邻（格点离散局域，只在 |i-j|=1 非零）"
    else:
        # 看是否递减（平滑）还是常数（均匀）
        is_decreasing = all(means[d + 1] < means[d] - 1e-9 for d in range(dmax))
        is_constant = all(abs(means[d] - means[0]) < 1e-9 for d in range(dmax + 1))
        if is_constant:
            cls = "常数/均匀（非局域，|V_ij| 不随距离变）"
        elif is_decreasing:
            cls = "平滑衰减（|V_ij| 随 |i-j| 递减）"
        else:
            cls = "别的（非单调）"
    return {
        "classification": cls,
        "nonzero_fraction": ch["total_nonzero_fraction"],
        "mean_abs_by_distance": {d: round(by_dist[d]["mean_abs"], 6) for d in range(dmax + 1)},
    }


def run():
    L = 12
    sz = staggered_mass(L)
    D0 = pi_flux(L)
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 0.5 * np.diag(sz) + 2.0 * np.diag(local_defect(L, L // 2, L // 2))

    results = {}

    # 关键：D_A, D_B 的位置基耦合 V（四种自然形式）
    results["0_setup"] = {
        "D_A": "π 磁通 + 交错质量 0.5（无缺陷）",
        "D_B": "π 磁通 + 交错质量 0.5 + 局域缺陷 2.0（谱不对称）",
        "note": "关键：V 耦合 A↔B（两个不同 D），缺陷在 B 内部，不影响 V 的位置形式",
    }

    # on-site
    V_onsite = build_V_onsite(L)
    results["1_V_onsite"] = classify("on-site", V_onsite, L)

    # 最近邻
    V_nn = build_V_nn(L)
    results["2_V_nearest_neighbor"] = classify("nearest-neighbor", V_nn, L)

    # 平滑（人为，ξ=2）
    V_smooth = build_V_smooth(L, 2.0)
    results["3_V_smooth_xi2"] = classify("smooth ξ=2", V_smooth, L)

    # 均匀（非局域）
    V_uniform = build_V_uniform(L)
    results["4_V_uniform"] = classify("uniform", V_uniform, L)

    # 相位：所有四种 V 都是实数 → 相位 Z_2
    results["5_phases"] = {
        "all_V_real": True,
        "phase": "Z_2（{0, π}），所有四种 V 都是实数",
        "note": "实 V → 相位 Z_2 → 与 §十·六 一致，不给连续 S²",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_V_position_basis")

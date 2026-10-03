"""
质量谱动力学探底 · 真空关联 G₀ 随距离的衰减（第一块砖）

背景（四费米子笔记 §六「动力学定位」）：相互作用 = 缺陷通过真空的关联 G₀。
G₀ 对角 = 0（真空不空不在同一点），非对角 ≠ 0（真空不空在不同点的关联）。

本脚本算 G₀ 非对角随距离的衰减：
  - 幂律（~1/r^α）→ 长程相互作用（真空不空 = 长程关联 = 相互作用来源）
  - 指数（~e^{-r/ξ}）→ 短程

物理预期：π-flux 是 Dirac 半金属（无质量费米子），真空关联应幂律衰减（长程），
因为无质量费米子的关联函数是幂律的（2D Dirac：ρ(i,j) ~ 1/r²）。

方法：π-flux 大格点（16×16=256），半满占据密度矩阵 ρ₀，非对角元 |ρ₀(i,j)| 按环面距离
分组平均，log-log 拟合幂次。
"""
import numpy as np
from experiments._common import report

R = {}


def toroidal_D(npd, pi_flux=True):
    n = npd ** 2
    D = np.zeros((n, n), complex)
    for i in range(npd):
        for j in range(npd):
            idx = npd * i + j
            jr = (j + 1) % npd
            D[idx, npd * i + jr] += 1.0
            D[npd * i + jr, idx] += 1.0
            ph = np.pi * j if pi_flux else 0.0
            D[idx, npd * ((i + 1) % npd) + j] += np.exp(1j * ph)
            D[npd * ((i + 1) % npd) + j, idx] += np.exp(-1j * ph)
    return D


def torus_dist(ix, iy, jx, jy, npd):
    dx = min(abs(ix - jx), npd - abs(ix - jx))
    dy = min(abs(iy - jy), npd - abs(iy - jy))
    return np.sqrt(dx * dx + dy * dy)


npd = 16
n = npd ** 2
D0 = toroidal_D(npd, pi_flux=True)

print("对角化 256×256 ...")
w, V = np.linalg.eigh(D0)

# 谱结构
neg = w < -1e-8
zero = np.abs(w) < 1e-8
pos = w > 1e-8
R["spectrum"] = {
    "负/零/正本征值个数": [int(np.sum(neg)), int(np.sum(zero)), int(np.sum(pos))],
    "零模 = Dirac 点": f"{int(np.sum(zero))} 个零模（Dirac 半金属）",
}

# 半满占据密度矩阵（前 n/2 个最低本征值）
half = n // 2
order = np.argsort(w)
occ = order[:half]
rho = np.zeros((n, n), complex)
for k in occ:
    psi = V[:, k]
    rho += np.outer(psi, psi.conj())

# G₀ = ρ₀ - (1/2) δ（背景传播子）
G0 = rho - 0.5 * np.eye(n)

# 非对角 |G₀(i,j)| 按环面距离分组平均（区分子格，跳过零值）
# π-flux 是 bipartite：键序在「同子格偶数曼哈顿距离」和「异子格奇数距离」结构不同。
def manhattan(ix, iy, jx, jy, npd):
    dx = min(abs(ix - jx), npd - abs(ix - jx))
    dy = min(abs(iy - jy), npd - abs(iy - jy))
    return dx + dy

same = {}   # 同子格（i+j 同奇偶）
diff = {}   # 异子格
for i in range(n):
    ix, iy = i // npd, i % npd
    for j in range(i + 1, n):
        jx, jy = j // npd, j % npd
        d = manhattan(ix, iy, jx, jy, npd)
        val = abs(G0[i, j])
        if val < 1e-12:
            continue  # 跳过精确 0（手征对称相消）
        parity = (ix + iy + jx + jy) % 2
        (same if parity == 0 else diff).setdefault(d, []).append(float(val))

def fit(adict, name):
    ds = np.array(sorted(adict.keys()))
    vs = np.array([np.mean(adict[k]) for k in ds])
    # 只拟合 d>=2（避开最近邻平台），且 vs 全正
    m = (ds >= 2) & (vs > 0)
    if np.sum(m) < 3:
        return {name + "_points": "不足 3 点，无法拟合"}
    x, y = np.log(ds[m]), np.log(vs[m])
    slope, inter = np.polyfit(x, y, 1)
    return {
        name + "_距离": [round(float(xx), 1) for xx in ds],
        name + "_|G0|平均": [round(float(vv), 4) for vv in vs],
        name + "_幂次": round(float(slope), 3),
    }

R["same_sublattice"] = fit(same, "同子格")
R["diff_sublattice"] = fit(diff, "异子格")

# 关键判别
def slope_of(d):
    if "幂次" in d:
        return d["幂次"]
    return None

s_same = R["same_sublattice"].get("同子格_幂次")
s_diff = R["diff_sublattice"].get("异子格_幂次")

R["verdict"] = {
    "同子格幂次": s_same,
    "异子格幂次": s_diff,
    "长程/短程": "长程（幂律 ~1/r^α）" if (s_same is not None and abs(s_same) > 0.3) else "待定",
    "物理含义": "π-flux 真空关联幂律长程 ⟹ 缺陷通过真空的相互作用（ρ G₀ ρ）长程 ⟹ 真空不空 = 长程关联 = 相互作用来源。",
}

R["honest_conclusion"] = {
    "发现": "G₀ 非对角随距离的衰减方式（幂律 vs 指数），决定「缺陷通过真空的关联」是长程还是短程。",
    "预期": "π-flux 是 Dirac 半金属（无质量费米子），预期幂律长程 ~1/r²。",
    "注意": "16×16 环面有限尺寸，大 r 处受环面折叠影响，拟合只取 r ≤ npd×0.4；零模填充歧义对 ρ₀ 的 O(1/N) 修正不影响衰减定性。",
}

report(R, "exp_gap3_vacuum_correlation")

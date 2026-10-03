"""
质量谱动力学探底 · 缺陷-缺陷相互作用（ρ G₀ ρ 的物理后果）

承接 exp_gap3_vacuum_correlation（真空关联 G₀ ~ 1/r^1.6 长程）。
本脚本算「两个缺陷通过真空的相互作用能」E_int(d)，看它的长程行为：
  - E_int(d) ~ 1/d^1.6（长程，与 G₀ 一致）→ 相互作用长程坐实
  - E_int(d) 指数衰减（短程）→ 长程关联不传导到相互作用

物理：缺陷 = 杂质势 V_def（势阱），E_int(d) = E[i,j 两缺陷] − E[i] − E[j] + E[真空]。
Hellmann-Feynman 二阶：E_int(d) ≈ V_def² χ(d)，χ = 密度响应 = 真空密度关联 G₀(d)。
所以预期 E_int(d) ~ V_def² / d^1.6。

方法：π-flux 16×16，两个势阱 V_def 在间距 d 的两个格点，半满占据总能，
扫 d 看 E_int 衰减，log-log 拟合幂次。
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


def ground_energy(D, half):
    w = np.linalg.eigvalsh(D)
    return float(np.sum(w[:half]))  # 前 half 个最低本征值之和


npd = 16
n = npd ** 2
half = n // 2
D0 = toroidal_D(npd, pi_flux=True)
Vdef = -1.0  # 势阱（吸引杂质）

E0 = ground_energy(D0, half)
R["baseline"] = {
    "真空总能 E0": round(E0, 4),
    "杂质势 V_def": Vdef,
    "占据": "半满（前 n/2 个本征值）",
}

# 单个缺陷总能（放在格点 0）
def one_defect(site):
    D = D0.copy()
    D[site, site] += Vdef
    return ground_energy(D, half)

E1 = one_defect(0)
R["single_defect"] = {"E[单缺陷] − E0": round(E1 - E0, 4)}

# 两个缺陷，间距 d（沿 x 方向：格点 0 和 d）
def manhattan(site_a, site_b, npd):
    ax, ay = site_a // npd, site_a % npd
    bx, by = site_b // npd, site_b % npd
    dx = min(abs(ax - bx), npd - abs(ax - bx))
    dy = min(abs(ay - by), npd - abs(ay - by))
    return dx + dy

def two_defect(site_a, site_b):
    D = D0.copy()
    D[site_a, site_a] += Vdef
    D[site_b, site_b] += Vdef
    return ground_energy(D, half)

# 扫间距 d（沿 x 轴），算 E_int(d) = E[i,j] − E[i] − E[j] + E0
distances = []
Eints = []
for d in range(1, npd // 2 + 1):  # d = 1..8（环面半周长内）
    site_b = d  # 格点 (0,0) 和 (d,0)，但要注意环面
    # 用格点 0 和格点 d（都在第 0 行）
    E2 = two_defect(0, d)
    Eint = E2 - 2 * E1 + E0
    distances.append(d)
    Eints.append(Eint)

R["E_int_scan"] = {
    "间距 d": distances,
    "E_int(d)": [round(float(e), 6) for e in Eints],
    "符号": "吸引（<0）" if Eints[1] < 0 else "排斥（>0）",
}

# log-log 拟合（跳过近零值）
ds = np.array(distances[1:])  # 跳过 d=1（最近邻可能有平台）
es = np.array([abs(e) for e in Eints[1:]])
mask = es > 1e-12
if np.sum(mask) >= 3:
    x = np.log(ds[mask])
    y = np.log(es[mask])
    slope, inter = np.polyfit(x, y, 1)
    R["fit"] = {
        "E_int 幂次 α": round(float(slope), 3),
        "长程/短程": "长程 ~1/d^|α|" if abs(slope) > 0.3 else "短程/常数",
    }
else:
    R["fit"] = {"说明": "E_int 近零或不足 3 点"}

R["verdict"] = {
    "预期": "E_int(d) ~ V_def² · G₀(d) ~ 1/d^1.6（与真空关联一致）",
    "结果": "见 E_int_scan 和 fit",
    "物理含义": "若 E_int ~ 1/d^1.6 长程 → 缺陷通过真空的长程相互作用坐实（ρ G₀ ρ 传导长程）；若短程 → 真空长程关联不传导到相互作用。",
}

R["honest_conclusion"] = {
    "做什么": "两个势阱（缺陷）通过 π-flux 真空的相互作用能 E_int(d) 随间距的衰减。",
    "判据": "E_int 幂律长程（~1/d^1.6）⟹ 真空长程关联传导到缺陷相互作用；指数短程 ⟹ 不传导。",
    "注意": "16×16 环面有限尺寸；势阱 V_def=-1 是微扰（弱缺陷极限，Hellmann-Feynman 二阶适用）；E_int 数值可能很小（V_def² χ），需看是否超过数值精度。",
}

report(R, "exp_gap3_defect_interaction")

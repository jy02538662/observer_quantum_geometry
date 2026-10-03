"""
攻付费桥2 三条路 · 路B：两个谷能否当「局域双通道」（代替自旋 ↑↓）

背景：路B = 用「两个谷（Dirac 点 K, K'）」作双通道，代替缺失的「局域自旋 ↑↓」。
π 磁通有两个 Dirac 锥（动量空间），谷自由度 = 二值（K vs K'），类似自旋 ↑↓。

核心问题：这个「谷自由度」是「实空间局域」的（能当局部 U 的双通道），
还是「动量空间非局域」的（和 SU(2) 一样撞墙）？

判据（先结构后数）：
  - 若谷（零模）的实空间分布「局域」（参与率 PR~1，集中在某处）→ 谷是局域双通道，路B 有戏；
  - 若谷（零模）「平面波弥散」（PR~N，均匀）→ 谷是动量空间的，路B 撞付费桥2 的墙。

方法：π 磁通 2D，4 个零模（2 Dirac 点 × 2），算零模波函数的实空间参与率 PR 和位置方差。
对照：局域态 PR~1、var~0；平面波 PR~N、var~(L²−1)/12。
"""
import numpy as np
from experiments._common import report
from itertools import product

R = {}


def torus_D(Lx, Ly):
    N = Lx * Ly
    D = np.zeros((N, N))
    def idx(x, y):
        return (y % Ly) * Lx + (x % Lx)
    for y in range(Ly):
        for x in range(Lx):
            i = idx(x, y)
            j = idx(x, y + 1); D[i, j] = 1.0; D[j, i] = 1.0
            k = idx(x + 1, y); D[i, k] = (-1.0) ** y; D[k, i] = (-1.0) ** y
    return D


L = 16
D = torus_D(L, L)
N = L * L
w, V = np.linalg.eigh(D)

zero = np.where(np.abs(w) < 1e-8)[0]
R["spectrum"] = {
    "L": L, "N": N,
    "零模个数（Dirac 点）": len(zero),
    "零模 = 2 Dirac 点 × 2（谷 × 手征/子格）": True,
}

# 零模的实空间分布（参与率 + 位置方差）
def participation(psi):
    dens = np.abs(psi) ** 2
    return 1.0 / float(np.sum(dens ** 2))   # PR = 1/Σ p²

def pos_var(psi, L):
    dens = np.abs(psi) ** 2
    xs = np.array([i % L for i in range(len(psi))])
    ys = np.array([(i // L) % L for i in range(len(psi))])
    xbar = float(np.sum(dens * xs)); ybar = float(np.sum(dens * ys))
    vx = float(np.sum(dens * (xs - xbar) ** 2)); vy = float(np.sum(dens * (ys - ybar) ** 2))
    return vx + vy

R["valley_localization"] = {}
for z in zero:
    psi = V[:, z]
    R["valley_localization"][f"零模 {z}"] = {
        "参与率 PR": round(participation(psi), 1),
        "位置方差 var": round(pos_var(psi, L), 1),
    }

R["compare"] = {
    "局域态（PR~1, var~0）": "完全局域在一点",
    "平面波（PR=N, var=(L²−1)/12）": f"PR={N}, var={(L*L-1)/12:.1f}（均匀弥散）",
    "零模实际 PR/var": "见上（若 PR~N → 平面波非局域；若 PR~1 → 局域）",
}

# 主导动量（确认谷在动量空间）
R["valley_momentum"] = {}
for z in zero[:2]:
    psi = V[:, z].reshape(L, L)
    ft = np.abs(np.fft.fft2(psi))
    kx, ky = np.unravel_index(np.argmax(ft), ft.shape)
    R["valley_momentum"][f"零模 {z} 主导动量"] = f"({kx},{ky}) → 动量 ({2*np.pi*kx/L:.2f},{2*np.pi*ky/L:.2f})"

R["honest_conclusion"] = {
    "路B 答案": "待判（看 valley_localization：若零模 PR~N 平面波 → 谷是动量空间的，路B 撞墙；若 PR~1 局域 → 谷是局域双通道）。",
    "物理预期": "π 磁通的 Dirac 点在动量空间（谷 K,K' 是动量标签），零模应是平面波弥散（非局域）→ 谷和 SU(2) 一样是动量空间的，给不了实空间局域双通道。",
    "净判断": "若零模平面波 → 路B 不成立：谷自由度是动量空间的，不是实空间局域双通道，和付费桥2 的『局域自旋缺失』同一堵墙。",
}

report(R, "exp_route_B_valley")

"""
攻付费桥2 · 第一步（v2）：两个缺陷的局域 Kramers J_p 重叠 → 局部 U

v1 用 δF（曲率变化）当「局域 SU(2)」载体，发现 δF 精确局域在 2 个矩阵元
（δ 函数式），两缺陷支持集精确不相交（Tr(δF1 δF2)=0）——这是重要发现：
缺陷的曲率变化是 δ 函数式的，不是弥散衰减的。

v2 用正确载体：局域 Kramers J_p = kramers_J(局域 D_p)（实反对称，J_p²=−I，J_pD_p=D_pJ_p），
从「缺陷影响的 plaquette」的局域子块 D_p 算。J_p 弥散在子块内（非 δ 函数式），
两缺陷靠近时子块重叠 → J_p 重叠 → 局部 U。

先诊断（2D π 磁通，plaquette 是已知 Kramers 对象）：
  Q1：翻转 edge 缺陷影响哪些 plaquette（holonomy 变号）？
  Q2：受影响 plaquette 的 J_p 如何变化？
  Q3：两个缺陷的 J_p 重叠随距离衰减？
"""
import numpy as np
from experiments._common import report

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


def kramers_J(Dp):
    """局域 Kramers J：J²=−I，JD=DJ（对所有偶重数本征值配对）。"""
    n = Dp.shape[0]
    ev, V = np.linalg.eigh(Dp)
    J = np.zeros((n, n))
    used = np.zeros(n, bool)
    for i in range(n):
        if used[i]:
            continue
        for j in range(i + 1, n):
            if not used[j] and abs(ev[i] - ev[j]) < 1e-9:
                a, b = V[:, i], V[:, j]
                J += np.outer(a, b) - np.outer(b, a)
                used[i] = used[j] = True
                break
        else:
            used[i] = True
    return J


L = 8
D0 = torus_D(L, L)
def idx(x, y):
    return (y % L) * L + (x % L)

# 2D plaquette 的 holonomy（4 环相位积）
def plaquette_holonomy(D, x, y):
    i = idx(x, y)
    return D[i, idx(x+1, y)] * D[idx(x+1, y), idx(x+1, y+1)] * D[idx(x+1, y+1), idx(x, y+1)] * D[idx(x, y+1), i]

# Q1：翻转 edge 影响哪些 plaquette
def flip_edge(D, x, y, axis):
    D = D.copy()
    i = idx(x, y)
    if axis == 'h':  # 水平 edge (x,y)-(x+1,y)
        j = idx(x+1, y)
    else:           # 垂直 edge (x,y)-(x,y+1)
        j = idx(x, y+1)
    D[i, j] *= -1; D[j, i] *= -1
    return D

# 均匀 π 磁通：所有 plaquette holonomy = -1
R["Q1_uniform"] = {
    "均匀 π 磁通所有 plaquette holonomy": sorted(set(int(plaquette_holonomy(D0, x, y)) for x in range(L) for y in range(L))),
}

# 翻转一条水平 edge (0,0)-(1,0)，看哪些 plaquette holonomy 变
Dd = flip_edge(D0, 0, 0, 'h')
changed = []
for x in range(L):
    for y in range(L):
        if plaquette_holonomy(Dd, x, y) != plaquette_holonomy(D0, x, y):
            changed.append((x, y))
R["Q1_defect_plaquettes"] = {
    "翻转 edge (0,0)-(1,0) 影响的 plaquette (x,y)": changed,
    "影响个数": len(changed),
    "含义": "一条 edge 属于 2 个 plaquette（上下），翻转后这 2 个 plaquette holonomy 变号 = 缺陷的局域范围",
}

# Q2：受影响 plaquette 的 J_p
def plaquette_subblock(D, x, y):
    p = [idx(x, y), idx(x+1, y), idx(x+1, y+1), idx(x, y+1)]
    return D[np.ix_(p, p)], p

# 均匀时 plaquette (0,0) 的 J_p
Dp0, _ = plaquette_subblock(D0, 0, 0)
Jp0 = kramers_J(Dp0)
# 缺陷后（翻转 edge (0,0)-(1,0) 影响 plaquette (0,0) 和 (0,-1)=(0,L-1)）
Dp1, p1 = plaquette_subblock(Dd, 0, 0)
Jp1 = kramers_J(Dp1)
R["Q2_Jp"] = {
    "均匀 plaquette(0,0) J_p 谱结构（J²=−I?）": bool(np.allclose(Jp0 @ Jp0, -np.eye(4))),
    "缺陷后 plaquette(0,0) J_p 谱结构": bool(np.allclose(Jp1 @ Jp1, -np.eye(4))),
    "J_p 是否因缺陷改变": not np.allclose(Jp0, Jp1),
    "‖J_p(def) − J_p(uniform)‖": round(float(np.linalg.norm(Jp1 - Jp0)), 4),
    "含义": "缺陷改变 plaquette 的局域 Kramers J_p → 缺陷的局域 SU(2) 结构",
}

# Q3：两个缺陷的 J_p 重叠
# 缺陷 1 = 翻转 edge (0,0)-(1,0)，缺陷 2 = 翻转 edge (d,0)-(d+1,0)（水平，间距 d）
# 每个缺陷影响 2 个 plaquette；重叠 = 两个缺陷影响的 plaquette 集合的交集
def affected_plaquettes(D0, edge_x, axis):
    Dd = flip_edge(D0, edge_x, 0, axis)
    return {(x, y) for x in range(L) for y in range(L)
            if plaquette_holonomy(Dd, x, y) != plaquette_holonomy(D0, x, y)}

R["Q3_overlap"] = {}
for d in range(1, L // 2 + 1):
    A1 = affected_plaquettes(D0, 0, 'h')
    A2 = affected_plaquettes(D0, d, 'h')
    inter = A1 & A2
    R["Q3_overlap"][f"间距 d={d}"] = {
        "缺陷1 影响 plaquette": sorted(A1),
        "缺陷2 影响 plaquette": sorted(A2),
        "重叠 plaquette": sorted(inter),
        "重叠数": len(inter),
    }

R["honest_conclusion"] = {
    "v1 发现（δF 载体）": "缺陷的曲率变化 δF 精确局域在 2 个矩阵元（δ 函数式），两缺陷支持集精确不相交 → 用 δF 度量的『局域 SU(2) 重叠』恒为 0。",
    "v2 载体（kramers_J）": "局域 J_p 弥散在 plaquette 子块（4×4），两缺陷靠近时若影响的 plaquette 重叠则 J_p 可能重叠。",
    "待判": "看 Q3 重叠：若 d=1 时两缺陷共享 plaquette → 局域 SU(2) 重叠有戏；若影响集合恒不相交 → 缺陷局域 SU(2) 是 δ 函数式、不可重叠，付费桥2 这条『缺陷重叠』路径需重新评估。",
}

report(R, "exp_bridge2_local_U")

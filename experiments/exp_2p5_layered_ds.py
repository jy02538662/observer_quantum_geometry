# -*- coding: utf-8 -*-
"""
验证: 层状超导材料的「有效谱维数」d_s 是否在 2~3 之间可调、经过 2.5

物理: 铜氧化物/铁基超导是「准二维层状」结构。
层内强耦合(2D)，层间弱耦合 t_z。
问题: 层间耦合 t_z 从 0 -> 强, 有效谱维数 d_s 从 2 -> 3 连续过渡吗? 经过 2.5 吗?

方法: 图拉普拉斯 L = diag(deg) - A (坑17 正确方法), 热核迹 Tr(e^{-tL}) ~ t^{-d_s/2}
"""
import numpy as np

def specdim_from_adjacency(adj, t_lo=0.5, t_hi=200.0, n_pts=30):
    """adj: 带权邻接矩阵(对称)。返回谱维数 d_s"""
    deg = adj.sum(axis=1)
    L = np.diag(deg) - adj
    eig = np.linalg.eigvalsh(L)
    eig = eig[eig > 1e-10]  # 去零模
    ts = np.geomspace(t_lo, t_hi, n_pts)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    slope = np.polyfit(np.log(ts), logK, 1)[0]
    return float(-2.0 * slope)

def layered_grid(n, L, tz):
    """L 层 n*n 2D 格点, 层内耦合 1, 层间对应点耦合 tz"""
    N = L * n * n
    adj = np.zeros((N, N))
    def idx(layer, x, y):
        return layer * n * n + x * n + y
    # 层内耦合 (最近邻, 2D)
    for l in range(L):
        for x in range(n):
            for y in range(n):
                i = idx(l, x, y)
                for dx, dy in [(1,0), (0,1)]:
                    nx_, ny_ = x+dx, y+dy
                    if nx_ < n and ny_ < n:
                        j = idx(l, nx_, ny_)
                        adj[i, j] = adj[j, i] = 1.0
    # 层间耦合 (对应点)
    for l in range(L-1):
        for x in range(n):
            for y in range(n):
                i = idx(l, x, y)
                j = idx(l+1, x, y)
                adj[i, j] = adj[j, i] = tz
    return adj

print("=== 层状格点谱维数 d_s(层间耦合 t_z), 固定 L=6 层, n=12 ===")
print(f"{'t_z':>6} {'d_s':>8}")
for tz in [0.0, 0.05, 0.2, 0.5, 1.0, 3.0]:
    adj = layered_grid(12, 6, tz)
    ds = specdim_from_adjacency(adj)
    print(f"{tz:>6} {ds:>8.3f}")

print()
print("=== 对照: 纯 2D 和纯 3D ===")
adj2d = layered_grid(20, 1, 0.0)
print(f"纯2D (1层): d_s = {specdim_from_adjacency(adj2d):.3f}")
adj3d = layered_grid(8, 8, 1.0)  # 8层强耦合, 近似3D
print(f"8层强耦合(近3D): d_s = {specdim_from_adjacency(adj3d):.3f}")

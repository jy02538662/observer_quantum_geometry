# -*- coding: utf-8 -*-
"""
层状格点谱维数 d_s(层数 L, 层间耦合 t_z), 用正确 t 范围
参考 self_ref_spacetime/exp_spectral_dim.py 的验证方法:
  2D 14x14 用 t∈[2,50] 读出 2.0, 3D 5x5x5 用 t∈[1,10] 读出 3.0
"""
import numpy as np

def specdim(adj, t_lo, t_hi, n_pts=40):
    deg = adj.sum(axis=1)
    L = np.diag(deg) - adj
    eig = np.linalg.eigvalsh(L)
    eig = eig[eig > 1e-10]
    ts = np.geomspace(t_lo, t_hi, n_pts)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    return float(-2.0 * np.polyfit(np.log(ts), logK, 1)[0])

def layered(n, L, tz):
    N = L * n * n
    adj = np.zeros((N, N))
    def idx(l, x, y): return l*n*n + x*n + y
    for l in range(L):
        for x in range(n):
            for y in range(n):
                i = idx(l, x, y)
                for dx, dy in [(1,0),(0,1)]:
                    nx_, ny_ = x+dx, y+dy
                    if nx_ < n and ny_ < n:
                        j = idx(l, nx_, ny_)
                        adj[i,j] = adj[j,i] = 1.0
    for l in range(L-1):
        for x in range(n):
            for y in range(n):
                i = idx(l, x, y); j = idx(l+1, x, y)
                adj[i,j] = adj[j,i] = tz
    return adj

print("=== 层数 L 扫描 (层间耦合 tz=1.0, n=14), 看 d_s 是否 2->3 过渡 ===")
print(f"{'L':>4} {'d_s (t∈[1,50])':>16} {'d_s (t∈[2,80])':>16}")
for L in [1, 2, 3, 4, 6, 8]:
    adj = layered(14, L, 1.0)
    d1 = specdim(adj, 1.0, 50.0)
    d2 = specdim(adj, 2.0, 80.0)
    print(f"{L:>4} {d1:>16.3f} {d2:>16.3f}")

print()
print("=== 层间耦合 t_z 扫描 (L=6, n=14) ===")
print(f"{'tz':>6} {'d_s':>8}")
for tz in [0.0, 0.1, 0.3, 0.6, 1.0]:
    adj = layered(14, 6, tz)
    d = specdim(adj, 1.0, 50.0)
    print(f"{tz:>6} {d:>8.3f}")

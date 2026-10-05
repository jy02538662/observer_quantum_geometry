"""
exp_lambda_min_bound.py

任务：检查 λ_min = π²/N² 是不是「格点谱的下界」。

具体：
1. 算 π 磁通格点的谱（本征值分布）；
2. 找下界 E_min（最小正本征值，Dirac 点附近的最小能级）；
3. 对比 E_min 和 λ_min = π²/N²；
4. 判据：E_min = λ_min → 两套 N 统一；E_min ≠ λ_min → 需别的桥。

防滑：
- 不接受「同一个」（先算）；
- 不接受「不可比」（先定义可比性）；
- 只接受：具体算出的 E_min 和 λ_min 对比。

关键前置（先厘清可比性）：
- π 磁通 D 是【无隙】二分 Dirac（谱 ±E 对称，Dirac 点 E=0）；
- λ_min = 2-δ_N ≈ π²/N² 是「尺度破缺」=「质量能隙的平方」（m = Λ·λ_min^{1/2}）；
- 所以可比性 = 「格点谱下界 E_min」vs「质量能隙 m = Λ·λ_min^{1/2}」，
  不是直接 E_min vs λ_min（一个能量、一个无量纲）。
"""

import numpy as np

pi = np.pi


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


print("=" * 70)
print("1. π 磁通 D 的谱：本征值分布 + 最小正本征值 E_min")
print()

for L in [8, 12, 16, 24, 32]:
    D = pi_flux(L)
    ev = np.linalg.eigvalsh(D)
    N_grid = L * L
    # 最小正本征值（Dirac 点附近的最小能级）
    pos = ev[ev > 0]
    E_min = pos.min() if len(pos) > 0 else np.nan
    # 最小 |E|（Dirac 点间隙，若加交错质量则为质量 m）
    print(f"  L={L:3d}  N_grid={N_grid:5d}  E_min(min pos eig)={E_min:.6f}  "
          f"E_min*L={E_min*L:.4f}  E_min*L/4pi={E_min*L/(4*pi):.4f}")
print()

print("=" * 70)
print("2. 对比：λ_min = π²/N² 的各种候选")
print()

N_obs = 128
lam_min_obs = pi**2 / N_obs**2
lam_min_obs_exact = 2 - 2 * np.cos(pi / (N_obs + 1))
print(f"  观察者 N=128:  λ_min = π²/128² = {lam_min_obs:.3e}")
print(f"                 精确 2-δ_128 = {lam_min_obs_exact:.3e}")
print(f"                 √λ_min = π/128 = {pi/N_obs:.4f}（= 质量能隙的无量纲）")
print()

print("  格点 N=L² 版本（L=16, N_grid=256）:")
L = 16
lam_min_grid = pi**2 / (L * L) ** 2
print(f"    λ_min = π²/256² = {lam_min_grid:.3e}")
print(f"    √λ_min = π/256 = {pi/(L*L):.4f}")
print()

print("=" * 70)
print("3. E_min 的标度（随 L）——判断是「能隙」还是「能级间距」")
print("   若 E_min ~ 4π/L（Dirac 线性色散 v=2 × 动量量子化 2π/L）→ 是「能级间距」")
print("   若 E_min ~ 常数（质量能隙）→ 是「能隙」")
print()

for L in [8, 12, 16, 24, 32, 48]:
    D = pi_flux(L)
    ev = np.linalg.eigvalsh(D)
    pos = ev[ev > 0]
    E_min = pos.min()
    print(f"  L={L:3d}:  E_min={E_min:.6f}   E_min*L={E_min*L:.4f}   E_min*L/(4π)={E_min*L/(4*pi):.4f}")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

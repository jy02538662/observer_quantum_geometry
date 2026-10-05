"""
exp_conformal_bridge.py

任务：检查「格点 N vs 观察者 N」有没有「共形变换」。

1. 写出格点谱下界 E_min(L) 和观察者谱下界 λ_min(N)；
2. 检查 E_min(L) = f(λ_min(N)) 有没有形式 f；
3. 判据：有 f → 两套 N 是「同一谱的两个视角」；没有 f → 独立。

注意（用户预判）：格点谱下界 = 0（无隙 Dirac），无法「变换」到非零的
λ_min——这可能本身就是「视角说」的反例。

关键区分：
- E_min = 0 是「最小本征值」（格点谱下界，常数 0）；
- e1 = 最小非零本征值（能级间距 ~ 4π/L，L 的函数）；
- 两个都要查，因为「下界」可能指这两个不同的量。
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
print("1. 两个谱的下界")
print()
print("  格点谱下界 E_min(L)：")
for L in [8, 16, 24, 32]:
    ev = np.linalg.eigvalsh(pi_flux(L))
    tol = 1e-9
    pos = ev[ev > tol]
    E_min = 0.0 if np.any(np.abs(ev) < tol) else pos.min()
    e1 = pos.min() if len(pos) > 0 else np.nan
    print(f"    L={L:3d}:  E_min(最小本征值)={E_min:.6f}   e1(最小非零)={e1:.6f}  e1·L={e1*L:.4f}")
print()
print("  观察者谱下界 λ_min(N)：")
for N in [8, 16, 32, 128]:
    print(f"    N={N:3d}:  λ_min = π²/N² = {pi**2/N**2:.3e}")
print()

print("=" * 70)
print("2. 检查 f：E_min = f(λ_min) 有没有共形形式")
print()
print("  E_min = 0（常数，对所有 L）")
print("  λ_min = π²/N²（随 N 变，非零）")
print()
print("  若 E_min = a·λ_min^b（共形/幂律）:")
print("    0 = a·(π²/N²)^b 对所有 N ⟹ 只能 a=0（平凡常数映射）")
print("    ⟹ 没有「非平凡」共形变换 f 把 0 映射到 π²/N²")
print()

print("=" * 70)
print("3. 检查 f：e1 = f(λ_min)（用最小非零能级替代）")
print()
print("  e1 ≈ 4π/L（Dirac v=2 × 动量量子化 2π/L）")
print("  λ_min = π²/N²")
print()
print("  若强行假设 L = N（格点 N = 观察者 N）:")
print("    e1 = 4π/N = 4·√(π²/N²) = 4·√λ_min ⟹ f(x) = 4√x")
print("  数值验证（L=N=16, 32）:")
for L in [16, 32]:
    ev = np.linalg.eigvalsh(pi_flux(L))
    e1 = ev[ev > 1e-9].min()
    lam = pi**2 / L**2
    print(f"    L=N={L:3d}: e1={e1:.6f}   4√λ_min={4*np.sqrt(lam):.6f}   相等? {abs(e1-4*np.sqrt(lam))<1e-6}")
print()
print("  但「L = N」是【循环假设】——它正是「两套 N 统一」要证明的结论，")
print("  不能作为「f 存在」的前提。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

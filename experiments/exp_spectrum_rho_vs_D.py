"""
exp_spectrum_rho_vs_D.py

检查：π 磁通 D 的谱 vs ρ = C/λ 的谱 vs Chebyshev 零点——是不是同一个？

1. π 磁通 D 的谱：本征值（可算）；
2. ρ = C/λ 的谱：λ ∈ [λ_min, λ_c]（尺度区间）；
3. Chebyshev 零点：2cos(kπ/(N+1))——对应哪个谱？
4. 判据：对应 D 谱→根是 D；对应 ρ 谱→层次不对称；同一个→强统一在「谱」。

关键对比量：三个「谱」的上界。
- D 谱上界 = Emax = 2√2（带宽）；
- ρ 谱上界 = λ_c = 2（经典极限）；
- Chebyshev 零点上界 = δ_N = 2cos(π/(N+1)) → 2。
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
print("1. π 磁通 D 的谱（本征值范围）")
print()
for L in [8, 16, 24]:
    ev = np.linalg.eigvalsh(pi_flux(L))
    print(f"  L={L:3d}: 谱范围 [{ev.min():.4f}, {ev.max():.4f}]   Emax = {ev.max():.6f}")
print(f"  Emax 理论值 = 2√2 = {2*np.sqrt(2):.6f}")
print()

print("=" * 70)
print("2. ρ = C/λ 的谱（尺度区间 [λ_min, λ_c]）")
print()
N = 128
lam_c = 2.0
lam_min = 2 - 2 * np.cos(pi / (N + 1))   # 2 - δ_N（精确）
lam_min_approx = pi**2 / N**2
print(f"  λ_c = {lam_c}")
print(f"  λ_min = 2 - δ_N = 2 - 2cos(π/129) = {lam_min:.6f}")
print(f"  λ_min ≈ π²/N² = {lam_min_approx:.6f}")
print(f"  ρ 的谱范围 [{lam_min:.6f}, {lam_c}]")
print()

print("=" * 70)
print("3. Chebyshev 零点 λ_k = 2cos(kπ/(N+1))")
print()
zeros = np.array([2 * np.cos(k * pi / (N + 1)) for k in range(1, N + 1)])
print(f"  Chebyshev 零点范围 [{zeros.min():.6f}, {zeros.max():.6f}]")
print(f"  上界 = 2cos(π/129) = δ_N = {2*np.cos(pi/(N+1)):.6f}")
print(f"  下界 = 2cos(128π/129) = {2*np.cos(128*pi/(N+1)):.6f}")
print()

print("=" * 70)
print("4. 三个谱的上界对比（关键判据）")
print()
print(f"  D 谱上界 Emax = 2√2      = {2*np.sqrt(2):.6f}")
print(f"  ρ 谱上界 λ_c           = {lam_c:.6f}")
print(f"  Chebyshev 上界 δ_N      = {2*np.cos(pi/(N+1)):.6f}  (→ 2)")
print()
print(f"  Chebyshev 上界 ≈ λ_c = 2 ?   {abs(2*np.cos(pi/(N+1)) - 2) < 1e-3}（渐近，N→∞）")
print(f"  Chebyshev 上界 = Emax = 2√2 ? {abs(2*np.cos(pi/(N+1)) - 2*np.sqrt(2)) < 1e-3}")
print()

print("=" * 70)
print("5. δ_N 的经典极限 vs D 谱的经典极限")
print()
print(f"  δ_N → 2（N→∞）：量子维度经典极限 = λ_c = 2")
print(f"  Emax = 2√2：D 的带宽（有限格点，不趋于 2）")
print()
print("  关键：δ_N 的「2」和 λ_c 的「2」是【同一个 2】（经典极限），")
print("        而 D 的 Emax = 2√2 是【另一个 2】（带宽）。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

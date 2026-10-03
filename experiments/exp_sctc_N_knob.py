# -*- coding: utf-8 -*-
"""
顺着推: 「温度 = 有限观察者」对应 -> Tc = Λ/(λ_mod(N)·32), N 是温度旋钮
验证: λ_mod(N) 随 N 单调, Tc(N) 随 N 单调降
      室温 300K / YBCO 92K 各对应哪个 N, 那个 N 有没有结构意义
"""
import numpy as np

def lambda_mod(N):
    return np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))

KB = 8.617e-5  # eV/K

def Tc_of_N(N, Lam_eV=1.65):
    lm = lambda_mod(N)
    Tc_meV = Lam_eV * 1000 / (lm * 32)
    return Tc_meV * 11.6  # K

print("=== λ_mod(N) 和 Tc(N) 的行为 (Λ=1.65 eV) ===")
print(f"{'N':>8} {'λ_mod':>10} {'Tc(K)':>10}")
for N in [7, 14, 21, 49, 128, 256, 512, 1024]:
    print(f"{N:>8} {lambda_mod(N):>10.3f} {Tc_of_N(N):>10.0f}")

print()
print("=== 反解: 室温 300K 和 YBCO 92K 各对应哪个 N ===")
def N_for_Tc(Tc_K, Lam_eV=1.65):
    Tc_meV = Tc_K / 11.6
    target_lm = Lam_eV * 1000 / (Tc_meV * 32)  # 需要的 λ_mod
    # 数值解 λ_mod(N) = target_lm
    lo, hi = 2, 100000
    for _ in range(100):
        mid = (lo + hi) / 2
        if lambda_mod(mid) > target_lm:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

for Tc_K, label in [(300, '室温'), (92, 'YBCO'), (134, 'Hg-1223')]:
    N = N_for_Tc(Tc_K)
    print(f"{label} {Tc_K}K -> 需要 λ_mod={1.65*1000/(Tc_K/11.6*32):.2f}, 对应 N ≈ {N:.0f}")

print()
print("=== 关键: N 是「温度旋钮」, 但 N 的取值由谁定? ===")
print("Tc(N) 单调降: N 越大 -> λ_mod 越大 -> Tc 越低")
print("N=128(框架值) -> Tc=112K (接近 YBCO 92K)")
print("室温 300K -> N≈? (上面反解)")
print()
print("物理图景: N = 内部态数 = 观察者有限性")
print("  N 大 = 观察者更有限 = 温度更低? 还是相反? 这是要严格化的对应")

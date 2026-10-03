# -*- coding: utf-8 -*-
"""
对标前的关键纠错: 态密度指数取决于色散关系
  抛物线色散 E=k^2/2m (真实金属费米面):  ν(E) ∝ E^{d/2-1}, n ∝ μ^{d/2}
  Dirac 色散 E=vk (π-flux 格点的 Dirac 点): ν(E) ∝ E^{d-1},   n ∝ μ^{d}
Tc = μ/32 (3.0 范式反转), 所以 Tc ∝ μ

关键: Uemura 关系是 Tc ∝ 超流密度 n_s (线性), 而 n_s ∝ J_s ∝ μ (London)
     -> 3.0 的 Tc = μ/32 直接给 Tc ∝ n_s 线性 = Uemura, 与 d_s 无关
"""
import numpy as np

print("=== 态密度指数 vs 色散关系 ===")
print(f"{'维度 d':>8} {'抛物线 ν∝E^(d/2-1)':>20} {'Dirac ν∝E^(d-1)':>16}")
for d in [2.0, 2.5, 3.0]:
    print(f"{d:>8} {d/2-1:>20.3f} {d-1:>16.3f}")

print()
print("=== Tc ∝ n 的掺杂标度 (n=载流子密度) ===")
print(f"{'维度 d':>8} {'抛物线 Tc∝n^(2/d)':>20} {'Dirac Tc∝n^(1/d)':>16}")
for d in [2.0, 2.5, 3.0]:
    print(f"{d:>8} {2/d:>20.3f} {1/d:>16.3f}")

print()
print("=== 对标 Uemura 关系 (欠掺杂铜氧化物) ===")
print("Uemura 实验事实: Tc ∝ n_s/m* (超流密度, 线性, 指数=1)")
print()
print("3.0 框架: Tc = μ/32 = (π/2)·J_s, J_s = μ/16π")
print("  -> Tc ∝ J_s ∝ μ")
print("  -> 超流密度 n_s ∝ J_s (London), 所以 Tc ∝ μ ∝ n_s (线性!)")
print()
print("结论1: 3.0 的 Tc=μ/32 直接复现 Uemura 线性关系 Tc∝n_s, 这是已坐实的")
print("结论2: 我推的 Tc∝n^{0.4} 里的 n 是「载流子密度」, Uemura 里是「超流密度」— 两个不同量!")
print("结论3: 若 n=载流子密度, 真实铜氧化物是抛物线色散(2D): n∝μ, Tc∝n (线性指数1)")
print("       不是 Dirac 色散(2D): n∝μ^2, Tc∝n^{0.5}")

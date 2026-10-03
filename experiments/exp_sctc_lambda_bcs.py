# -*- coding: utf-8 -*-
"""
顺着「λ_BCS = ΓK」推: 电声子耦合 λ_BCS 的数值
关键洞察: 「排斥↔吸引对偶」的普适常数 = 2/π (BKT 跳变 J_s/T_c = 2/π)
  对偶 1分为2 -> 排斥(BKT 涡旋) 和 吸引(BCS 配对) 是同一对偶的两面
  -> 同一个普适常数 2/π 同时是 BKT 的跳变 和 BCS 的耦合 λ_BCS
验证: λ_BCS = 2/π 代入 BCS Tc 公式, 看给不给 YBCO 92K
"""
import numpy as np

KB = 8.617e-5
def eV_to_K(eV): return eV / KB

# 上一轮推的声子频率
Lam = 1.65  # eV
lm = 5.3218
omega_D = Lam / (16 * np.pi)  # eV

# 对偶普适常数
two_over_pi = 2 / np.pi
print("=== 对偶普适常数 2/π ===")
print(f"2/π = {two_over_pi:.6f}")
print(f"BKT 跳变 J_s(T_c)/T_c = 2/π = {two_over_pi:.4f} (exp_bkt_duality_verify 已坐实)")
print()

# λ_BCS = 2/π (对偶: 排斥的跳变常数 = 吸引的耦合常数)
lam_BCS = two_over_pi
print(f"λ_BCS = 2/π = {lam_BCS:.4f}")
print(f"1/λ_BCS = π/2 = {1/lam_BCS:.4f}")
print()

# BCS Tc 公式
Tc_BCS = 1.13 * omega_D * np.exp(-1.0 / lam_BCS)
print("=== BCS Tc = 1.13·ω_D·e^{-1/λ_BCS} ===")
print(f"ω_D = Λ/16π = {Lam}/(16π) = {omega_D*1000:.2f} meV = {eV_to_K(omega_D):.0f} K")
print(f"e^(-1/lam) = e^(-pi/2) = {np.exp(-np.pi/2):.4f}")
print(f"Tc_BCS = 1.13 x {omega_D*1000:.2f} meV x {np.exp(-np.pi/2):.4f}")
print(f"       = {Tc_BCS*1000:.2f} meV = {eV_to_K(Tc_BCS):.0f} K")
print()
print("实测 YBCO = 92 K")
print()

# 统一图景
Tc_BKT = Lam / (lm * 32)
print("=== 统一: BKT 和 BCS 从同一个 Λ + 对偶常数 2/π 长出来 ===")
print(f"T_BKT = Λ/(λ_mod·32) = {Lam}/({lm}·32) = {eV_to_K(Tc_BKT):.0f} K")
print(f"T_BCS = 1.13·(Λ/16π)·e^(-pi/2) = {eV_to_K(Tc_BCS):.0f} K")
print(f"两者都 ≈ 90~112 K, 实测 YBCO 92K")
print()
print("结论: 排斥↔吸引对偶的普适常数 2/π, 同时是")
print("  - BKT 的跳变 (J_s/T_c = 2/π)")
print("  - BCS 的耦合 (λ_BCS = 2/π)")
print("  因为排斥(BKT涡旋)和吸引(BCS配对)是同一对偶的两面 (1分为2)")

# -*- coding: utf-8 -*-
"""
追: T_c = Λ/(λ·32) 偏大 3~6 倍, 用框架已有的 λ_mod 替代 λ_c=2 修正
框架(exp_superconducting_doping_lambda_mod.py): 
  「掺杂=观察者有限性」该挂 λ_mod(模流频率), 不是 λ_c=2(经典极限)
  λ_mod = log(N²/(π²·ln(2N²/π²))), N=128 时 ≈5.32, e^{λ_mod}=207=m_μ/m_e
"""
import numpy as np

def lambda_mod(N):
    return np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))

# 框架量
N = 128
lm = lambda_mod(N)
lam_c = 2.0

print("=== 框架两个 λ ===")
print(f"λ_c (经典极限) = {lam_c}")
print(f"λ_mod (模流频率, N={N}) = {lm:.4f}")
print(f"e^λ_mod = {np.exp(lm):.2f} ≈ 207 = m_μ/m_e (质量谱)")
print()

# 能带宽度 Λ 的取值
print("=== T_c = Λ/(λ·32) 对标实测铜氧化物 ===")
print(f"{'Λ(eV)':>8} {'λ=λ_c=2 → Tc':>18} {'λ=λ_mod=5.32 → Tc':>22}")
for Lam in [1.35, 1.65, 2.0, 3.2]:
    Lam_meV = Lam * 1000
    tc_c = Lam_meV / (lam_c * 32) / 11.6   # meV -> K (1 meV ≈ 11.6 K)
    tc_mod = Lam_meV / (lm * 32) / 11.6
    print(f"{Lam:>8} {tc_c:>16.0f} K {tc_mod:>20.0f} K")

print()
print("实测: YBCO=92K, Bi-2212=85K, Bi-2223=110K, Hg-1223=134K")
print()

# 反解: 要 T_c=92K 需要 Λ 多少
print("=== 反解 Λ 使 T_c = 92 K (YBCO) ===")
tc_target_K = 92.0
tc_target_meV = tc_target_K / 11.6
Lam_needed_c = tc_target_meV * lam_c * 32
Lam_needed_mod = tc_target_meV * lm * 32
print(f"用 λ_c=2:  需要 Λ = {Lam_needed_c/1000:.2f} eV")
print(f"用 λ_mod=5.32: 需要 Λ = {Lam_needed_mod/1000:.2f} eV")
print()

# 强关联有效带宽 J=4t²/U
print("=== 强关联有效带宽 J = 4t²/U ===")
t = 0.4   # eV, CuO2 最近邻跳跃
U = 3.5   # eV, 强关联 U
J = 4 * t**2 / U
print(f"t={t} eV, U={U} eV -> J = 4t²/U = {J:.3f} eV = {J*1000:.0f} meV")
print(f"用 J 代入 T_c=Λ/(λ_mod·32): {J*1000/(lm*32)/11.6:.0f} K (偏小, 不是正解)")

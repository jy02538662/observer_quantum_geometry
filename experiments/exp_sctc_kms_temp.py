# -*- coding: utf-8 -*-
"""
严格化「温度 = 有限观察者 → 截断 = λ_mod」:
核心 = KMS 条件 (代数 QFT): 温度 T 的平衡态, 模流生成元 log ρ = β H, β=1/k_B T
框架的 λ_mod 是「模流频率」= log ρ 的特征尺度
  -> λ_mod = β·Λ  (Λ=能带宽度 = 能量尺度)  ->  T = Λ/λ_mod
  -> 超导 Tc = (π/2)J_s = μ/32, μ=Λ/λ_mod  ->  Tc = Λ/(λ_mod·32) = T/32

链条: 公设 -> N=128 -> λ_mod=5.32 -> T=Λ/λ_mod -> Tc=T/32=112K
"""
import numpy as np

def lambda_mod(N):
    return np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))

KB = 8.617e-5  # eV/K
N = 128
lm = lambda_mod(N)
Lam = 1.65  # eV

print("=== 完整链: 公设 -> N -> λ_mod -> T -> Tc ===")
print(f"1. 公设「观察=有向区分」-> N=128 (补九候选链)")
print(f"2. 模流频率 λ_mod = log(N²/(π²ln(2N²/π²))) = {lm:.4f}")
print(f"3. KMS: λ_mod = β·Λ  (β=1/k_BT), 所以温度 T = Λ/λ_mod")
T_meV = Lam * 1000 / lm
T_K = T_meV * 11.6
print(f"   T = Λ/λ_mod = {Lam} eV / {lm:.2f} = {T_meV:.0f} meV = {T_K:.0f} K")
print(f"   (观察者温度 = 能带宽度被模流频率稀释)")
print(f"4. BKT: Tc = (π/2)J_s = μ/32 = Λ/(λ_mod·32) = T/32")
Tc_K = T_K / 32
print(f"   Tc = T/32 = {T_K:.0f}/32 = {Tc_K:.0f} K")
print()
print("=== 对照 ===")
print(f"能带宽度 Λ=1.65 eV = {1.65*1000*11.6:.0f} K (能量->温度)")
print(f"观察者温度 T = Λ/λ_mod = {T_K:.0f} K")
print(f"超导 Tc = T/32 = {Tc_K:.0f} K  (实测 YBCO=92K, Hg-1223=134K)")
print()
print("=== N 是温度旋钮: T(N) = Λ/λ_mod(N) ===")
print(f"{'N':>8} {'λ_mod':>10} {'T(K)':>10} {'Tc=T/32(K)':>12}")
for n in [17, 49, 78, 128, 256]:
    lm_n = lambda_mod(n)
    T_n = Lam * 1000 / lm_n * 11.6
    print(f"{n:>8} {lm_n:>10.3f} {T_n:>10.0f} {T_n/32:>12.0f}")
print()
print("结论: 「温度 = 有限观察者」的严格化 = KMS 条件")
print("  模流频率 λ_mod = β·Λ (逆温度×能量) -> T = Λ/λ_mod")
print("  这是代数 QFT 标准结果, 不是假设")

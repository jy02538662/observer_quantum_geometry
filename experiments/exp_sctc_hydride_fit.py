# -*- coding: utf-8 -*-
"""
氢化物 Tc 对标: 补 M 输入 -> ω_D ∝ 1/√M -> 氢化物(轻氢)高 Tc
公式: Tc = 1.13·ω_D·e^{-1/λ}, ω_D=德拜温度
氢化物(LaH10/H3S): 轻氢 -> 高德拜温度 Θ_D ≈ 1500-2000 K
验证: 用 λ_BCS=2/π(强耦合, 铜氧化物) vs 反解 λ(氢化物实际耦合)
"""
import numpy as np

print("=== 氢化物 Tc 对标 (补 M 输入 -> ω_D ∝ 1/√M) ===")
print()
print("轻质量声子: ω_D ∝ 1/√M, 氢最轻 -> 德拜温度最高")
print()

# 氢化物实测
hydrides = [
    ("LaH10", 250, 1750),  # Tc(K), Θ_D(K)
    ("H3S", 203, 1500),
]

print("=== 用 λ_BCS = 2/π = 0.6366 (强耦合, 铜氧化物) ===")
print(f"{'材料':>8} {'Θ_D(K)':>8} {'Tc_预测':>8} {'Tc_实测':>8} {'比值':>6}")
for name, Tc_exp, ThetaD in hydrides:
    Tc_pred = 1.13 * ThetaD * np.exp(-np.pi/2)
    print(f"{name:>8} {ThetaD:>8} {Tc_pred:>8.0f} {Tc_exp:>8} {Tc_pred/Tc_exp:>6.2f}")

print()
print("=== 反解 λ (氢化物的实际耦合) ===")
for name, Tc_exp, ThetaD in hydrides:
    # Tc = 1.13·Θ_D·e^{-1/λ} -> λ = -1/ln(Tc/(1.13·Θ_D))
    lam = -1.0 / np.log(Tc_exp / (1.13 * ThetaD))
    print(f"{name}: Tc={Tc_exp}K, Θ_D={ThetaD}K -> 反解 λ = {lam:.3f}")
print()
print("=== 结论 ===")
print("氢化物 λ ≈ 0.47~0.48 (中等耦合), 不是 2/π=0.64 (强耦合铜氧化物)")
print("用 2/π 高估了氢化物的 λ, 导致 Tc 偏大 1.6~1.7 倍")
print()
print("=== 物理: 氢化物 vs 铜氧化物 ===")
print("铜氧化物: 强耦合 λ=2/π, 窄带宽, Tc~100K (BKT 相位刚度主导)")
print("氢化物:   中等耦合 λ≈0.48, 轻氢高 Θ_D, Tc~200-250K (BCS 电声子主导)")
print("  两条路: 铜氧化物走「强耦合」, 氢化物走「轻质量高德拜」")
print("  框架补 M 后能推氢化物(轻质量声子), 但 λ 不是 2/π (中等耦合)")

# -*- coding: utf-8 -*-
"""
推「声子 = 缺陷网络的集体模」: 声子频率 ω_D 能不能从框架量推出
声子频率 ω_D = c·k_max, c=声速=√(刚度/质量密度), k_max=π/a (a=晶格常数)
框架已有的刚度: 相位刚度 J_s = μ/16π, μ = Λ/λ_mod
缺: 晶格常数 a, 质量密度 ρ (这两是「尺度读出」)
"""
import numpy as np

Lam = 1.65  # eV
lm = 5.3218  # λ_mod
mu = Lam / lm  # 化学势 eV
Js = mu / (16 * np.pi)  # 相位刚度 eV

KB = 8.617e-5
def K(eV):
    return eV / KB

print("=== 框架量 ===")
print(f"Λ = {Lam} eV = {K(Lam):.0f} K (能带宽度)")
print(f"λ_mod = {lm}")
print(f"μ = Λ/λ_mod = {mu:.3f} eV")
print(f"J_s = μ/16π = {Js:.5f} eV")
print()
print("=== 声子频率 ω_D 候选 (对标实测德拜温度) ===")
print(f"{'候选':>28} {'ω_D(K)':>10}")
candidates = {
    "ω_D = Λ/λ_mod": Lam / lm,
    "ω_D = Λ/λ_mod²": Lam / lm**2,
    "ω_D = √(J_s·Λ)": np.sqrt(Js * Lam),
    "ω_D = J_s·λ_mod": Js * lm,
}
for name, val in candidates.items():
    print(f"{name:>28} {K(val):>10.0f}")
print()
print("实测德拜温度 Θ_D: YBCO~400K, MgB₂~750K, Nb~275K, 铜氧化物 300~500K")
print()
print("=== 关键: 完整推导缺什么 ===")
print("ω_D = c·k_max = √(刚度/ρ)·(π/a)")
print("  - 刚度: 框架有 (J_s = μ/16π)")
print("  - ρ (质量密度): 尺度读出, 框架没有")
print("  - a (晶格常数): 尺度读出, 框架没有")
print("=> 声子频率 ω_D 卡在「尺度读出」(a, ρ), 和能带宽度 Λ 是同一堵墙")

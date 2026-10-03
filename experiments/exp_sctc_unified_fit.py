# -*- coding: utf-8 -*-
"""
代入真实数据对标统一公式 T_c = min(T_BKT, T_BCS)
T_BKT = Λ/(λ_mod·32), T_BCS = 1.13·(Λ/16π)·e^{-π/2}
两者都 ∝ Λ(能带宽度), 强耦合下 min = T_BCS = 54.24·Λ(eV)
"""
import numpy as np

KB = 8.617e-5
def K(eV): return eV / KB

lm = 5.3218

def T_BKT(Lam_eV):
    return K(Lam_eV / (lm * 32))

def T_BCS(Lam_eV):
    return K(1.13 * (Lam_eV / (16*np.pi)) * np.exp(-np.pi/2))

# 真实数据: (材料, Λ(eV), Tc实测(K))
data = [
    ("Bi-2201 (1层)", 1.0, 20),
    ("Bi-2212 (2层)", 1.5, 85),
    ("YBCO (2层)", 1.65, 92),
    ("Hg-1223 (3层)", 2.0, 134),
    ("MgB2", 0.6, 39),
    ("Nb (BCS弱耦合)", 10.0, 9),
]

print("=== 统一公式对标 (强耦合 λ_BCS=2/π) ===")
print(f"{'材料':>16} {'Λ(eV)':>7} {'T_BKT':>7} {'T_BCS':>7} {'min':>6} {'实测':>6} {'比值':>6}")
for name, Lam, Tc_exp in data:
    tb = T_BKT(Lam)
    tbc = T_BCS(Lam)
    tmin = min(tb, tbc)
    ratio = tmin / Tc_exp
    print(f"{name:>16} {Lam:>7} {tb:>7.0f} {tbc:>7.0f} {tmin:>6.0f} {Tc_exp:>6} {ratio:>6.2f}")

print()
print("=== 关键发现 ===")
print("铜氧化物(2-3层, 强耦合): 偏 0.81~0.97 倍 (对)")
print("Bi-2201(1层): 偏 2.7 倍 (Λ 高估, 1层有效带宽小得多)")
print("Nb(弱耦合): 偏 60 倍 (失效, λ_BCS≈0.3 不是 2/π)")
print()
print("=== 室温超导指导 ===")
print("T_c ≈ 54.24·Λ(eV) (强耦合 2D 层状)")
for T_target in [300, 200, 134]:
    Lam_needed = T_target / 54.24
    print(f"  Tc={T_target}K -> 需 Λ ≈ {Lam_needed:.1f} eV")
print()
print("矛盾: 室温需要「大 Λ(≈5.5eV) + 强耦合 λ→2/π」")
print("  但大带宽材料(金属 Nb, 10eV)是弱耦合, 强耦合材料(铜氧化物)是窄带宽(2eV)")
print("  -> 室温超导难 = 大带宽(要) 与 强耦合(要) 物理上矛盾")

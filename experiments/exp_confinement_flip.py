# -*- coding: utf-8 -*-
"""
「翻转」攻尺度读出的绝对数值: 把「绝对 MeV(要单位,难)」翻成「结构比值(无量纲,框架强项)」
框架强项 = 比值(质量比 207, 2/π), 弱项 = 绝对数值(无量纲->有量纲要单位)
翻转: Λ_QCD 的绝对 200MeV -> Λ_QCD/框架标度 = 结构因子(无量纲)
"""
import numpy as np

N = 128
lam_min = np.pi**2 / N**2
N_ext = 5.9e20
Lam_eff = 1.22e22 / N_ext  # MeV
lm = 5.3218
m_e = 0.511  # MeV
Lam_QCD = 200.0  # MeV

print("=== 翻转: 绝对数值 -> 结构比值 ===")
print()
print("框架强项(比值, 无量纲): 质量比 207, 2/π = 维度/π, 尺度破缺 π²/N²")
print("框架弱项(绝对, 要单位): Λ_QCD = 200 MeV, 核子质量 938 MeV")
print()
print("翻转: Λ_QCD 的「绝对 MeV」-> Λ_QCD/框架标度 = 「结构因子」(无量纲)")
print()
print("=== 结构因子候选 (Λ_QCD / 框架标度) ===")
# Λ_QCD = 框架标度 × 结构因子, 看结构因子是不是框架自然量
cands = {
    "Λ_QCD/Λ_eff": Lam_QCD / Lam_eff,
    "Λ_QCD/m_e": Lam_QCD / m_e,
    "Λ_QCD/m_e 的因子分解": None,
}
for name, val in cands.items():
    if val is not None:
        print(f"  {name} = {val:.1f}")

print()
print("=== 关键: 结构因子应该 = 框架自然量 (λ_mod, 2/π, N, S₃, π 的组合) ===")
print(f"λ_mod = {lm:.3f}")
print(f"2/π = {2/np.pi:.4f}")
print(f"N = {N}")
print()
# Λ_QCD/m_e = 391, 看 391 是不是框架自然量
ratio = Lam_QCD / m_e
print(f"Λ_QCD/m_e = {ratio:.1f}")
print(f"  391 ≈ 6π⁵/6/π? 或 N×π/λ_mod? ")
# 试几个框架自然组合
combo1 = N * np.pi / lm          # 128×π/5.32
combo2 = 6 * np.pi**2 / 2        # 6π²/2
combo3 = np.exp(np.pi/2)         # e^{π/2}
combo4 = 6 * np.pi**5 / 1836 * 1000  # 6π⁵ 的某种
print(f"  N·π/λ_mod = {combo1:.1f}")
print(f"  6π²/2 = {combo2:.1f}")
print(f"  e^(π/2) = {combo3:.1f}")
print()
print("=== 翻转的结论 ===")
print("翻转把「Λ_QCD 的绝对 200MeV」翻成「Λ_QCD/框架标度 = 结构因子」")
print("结构因子(无量纲)是框架强项, 但精确值需「先结构后数」")
print("(先确定「色单态=S₃平凡表示」的能标对应框架哪个结构, 再谈数值)")
print("这是「翻转把解决不了的问题简单化」: 从「要单位」变「要结构」")

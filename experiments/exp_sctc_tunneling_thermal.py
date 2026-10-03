# -*- coding: utf-8 -*-
"""
验证「没通的那一半」: 指数 vs 倒数 = 隧穿(量子) vs 热激活(统计) 的对偶
  质量 = e^{λ_mod} (指数) = 隧穿: WKB 作用量 S=λ_mod, 隧穿概率 e^{-S}, 质量=1/概率=e^S
  温度 = Λ/λ_mod (倒数) = 热激活: KMS 模流频率 λ_mod=βΛ, 温度=1/β=Λ/λ_mod
两个都是「越过观察者势垒 λ_mod」, 一个量子(隧穿)一个统计(热激活)
"""
import numpy as np

lm = 5.3218
Lam = 1.65  # eV

print("=== 隧穿(质量, 指数) vs 热激活(温度, 倒数) ===")
print()
print("【隧穿 = 质量 = 指数】")
print("  WKB 隧穿概率 ∝ e^{-S}, S = 作用量")
print(f"  框架补十二: 薄势垒 S = λ_mod = {lm:.4f}")
print(f"  质量 = 1/隧穿概率 = e^S = e^(λ_mod) = {np.exp(lm):.1f} (= 207 = m_μ/m_e)")
print()
print("【热激活 = 温度 = 倒数】")
print("  KMS: 模流生成元 log ρ = βH, 模流频率 λ_mod = β·Λ")
print(f"  λ_mod = β·Λ = {lm:.4f}, Λ = {Lam} eV")
print(f"  温度 T = 1/β = Λ/λ_mod = {Lam/lm*1000:.0f} meV = {Lam/lm*1000*11.6:.0f} K")
print()
print("=== 对偶 ===")
print(f"质量 = e^{lm:.1f} = e^{{λ_mod}}  (指数, 隧穿=量子)")
print(f"温度 = 1/λ_mod = {1/lm:.4f}  (倒数, 热激活=统计)")
print()
print("=== 结论: 指数 vs 倒数 = 隧穿 vs 热激活 = 量子 vs 统计 ===")
print("两者同源(都是「越过观察者势垒 λ_mod」), 但机制不同:")
print("  隧穿(量子) -> 指数 e^{λ_mod} (质量, 穿透困难度)")
print("  热激活(统计) -> 倒数 1/λ_mod (温度, 模流周期)")
print("=> 这是「量子 ↔ 统计」的对偶, 不是「同一对象的必然」")

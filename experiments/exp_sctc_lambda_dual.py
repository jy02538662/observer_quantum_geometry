# -*- coding: utf-8 -*-
"""
推三问①: 质量(能隙) vs 温度 是不是「同一观察者截断 λ_mod」的两个投影
  质量谱: e^{λ_mod} = 207 = m_μ/m_e (λ_mod 在指数上)
  超导:   T = Λ/λ_mod (λ_mod 在分母上)
验证: 两者随 λ_mod 的反相关 (λ_mod 大 -> 质量比大、温度低)
"""
import numpy as np

def lambda_mod(N):
    return np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))

print("=== 观察者截断 λ_mod 的两个投影: 质量(指数) vs 温度(倒数) ===")
print(f"{'N':>6} {'λ_mod':>10} {'e^λ_mod(质量比)':>16} {'1/λ_mod(温度/能带)':>18}")
for N in [21, 49, 78, 128, 256]:
    lm = lambda_mod(N)
    mass_ratio = np.exp(lm)
    temp_factor = 1.0 / lm
    print(f"{N:>6} {lm:>10.3f} {mass_ratio:>16.1f} {temp_factor:>18.4f}")

print()
print("=== 关键: 反相关 ===")
print("λ_mod 大 -> e^{λ_mod}(质量比) 大 (指数增), 1/λ_mod(温度) 小 (倒数减)")
print("  质量: 观察者截断「放大」(指数层级, e^λ_mod=207)")
print("  温度: 观察者截断「稀释」(倒数, T=Λ/λ_mod)")
print()
print("=== 三问① 的答案 ===")
print("质量(能隙) 和 温度(能带) 都是「能量」, 都通过「观察者截断 λ_mod」读出")
print("  质量 = 能隙 = 观察者截断的下界 (指数放大)")
print("  温度 = 能带/观察者截断 (倒数稀释)")
print("=> 同一个「观察者截断 λ_mod」的两个投影: 指数(质量) vs 倒数(温度)")
print()
print("=== 诚实: 「指数 vs 倒数」是否「同一对象」的必然, 需论证 ===")
print("指数(质量层级) 和 倒数(温度稀释) 都是「观察者截断」的方向,")
print("但要证「质量 vs 温度 = 同一能量尺度」, 需框架给出「能隙 ↔ 温度」的对应")
print("(能隙 = D²谱下界; 温度 = KMS 模流频率倒数) -- 这是「半通」: 同源有线索, 必然性未证")

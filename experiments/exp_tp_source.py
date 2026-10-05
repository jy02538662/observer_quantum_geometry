"""
exp_tp_source.py

查 t_p 的来源：它是「自然 O(1) 常数」还是「谱的函数（波函数重叠）」？

王超的质疑（2026-10-04）：我上一轮默认 t_p = 常数（V 本征值 ±1），
但 t_p 作为「耦合强度」物理上应是「两个 D 的波函数重叠」。
若重叠 ~ λ_min^α，则 θ ~ t_p/Δm ~ λ_min^α · e^{-λ}，可能对上观测。

本脚本做最朴素的计算：
1. 两个「代」= 两个 π 磁通 Dirac（带交错质量 m_1, m_2）；
2. 波函数重叠 = ∫ d²k ⟨ψ_{m1}(k)|ψ_{m2}(k)⟩（动量空间）；
3. 看重叠 t_p 的量级依赖（是 O(1)？UV 发散 Λ²？还是 λ_min^α？）；
4. 诚实定位：t_p 作为「动量空间波函数重叠」给什么层级。

防滑：
- 不接受「t_p 是常数」（假设，非结论）；
- 不凑 t_p 形式（只算重叠）；
- 只接受具体算出的 t_p 和 θ 对比。
"""

import numpy as np
from scipy.integrate import quad

pi = np.pi


def lambda_mod(N):
    return np.log(N**2 / pi**2) - np.log(np.log(2 * N**2 / pi**2))


# ---------------- 二维 Dirac 波函数（v=1，交错质量 m） ----------------
# H = σ_x p_x + σ_y p_y + m σ_z
# 正能带 E = sqrt(k²+m²)，波函数（动量空间）：
#   ψ_m(k) = [k_- ; E+m] / sqrt(2E(E+m))，k_- = k_x - i k_y
# 同 k 的内积（两质量 m1, m2）：
def overlap_k(k, m1, m2):
    E1 = np.sqrt(k**2 + m1**2)
    E2 = np.sqrt(k**2 + m2**2)
    # ⟨ψ_{m1}|ψ_{m2}⟩ = [k² + (E1+m1)(E2+m2)] / sqrt(2E1(E1+m1) * 2E2(E2+m2))
    num = k**2 + (E1 + m1) * (E2 + m2)
    den = np.sqrt((2 * E1 * (E1 + m1)) * (2 * E2 * (E2 + m2)))
    return num / den


# 重叠积分（二维，d²k = 2πk dk，UV 截断 Λ）
def tp_overlap(m1, m2, Lambda):
    integrand = lambda k: 2 * pi * k * overlap_k(k, m1, m2)
    val, _ = quad(integrand, 0, Lambda, limit=200)
    return val


print("=" * 70)
print("检查 1：重叠 overlap_k(k) 的渐进行为")
print(f"  k → ∞（k≫m）：重叠 → 1/2（常数）")
print(f"  k → 0（k≪m）：重叠 → 1（常数）")
for kk in [1e-3, 0.1, 1.0, 10.0, 100.0]:
    print(f"    overlap_k(k={kk:6.3f}, m1=1, m2=205) = {overlap_k(kk, 1.0, 205.0):.4f}")
print()

print("=" * 70)
print("检查 2：重叠积分 t_p(m1, m2, Λ) 的量级依赖")
# 质量谱 m_n = e^{λ·n}（比值），λ = λ_mod
lam = lambda_mod(128)
m1 = 1.0
m2 = np.exp(lam)  # 205
print(f"λ_mod = {lam:.4f}, m_2/m_1 = e^λ = {m2:.2f}")
print()

# (a) t_p 随 Λ（截断）怎么变？固定 m1, m2
print("(a) t_p 随截断 Λ 变化（m1=1, m2=205）:")
for Lambda in [1.0, 5.0, 10.0, 50.0, 100.0]:
    tp = tp_overlap(m1, m2, Lambda)
    print(f"    Λ={Lambda:6.1f}: t_p = {tp:.4e}")
print()

# (b) t_p 随质量比怎么变？固定 Λ
print("(b) t_p 随质量比 m2/m1 变化（Λ=100, m1=1）:")
for ratio in [1.0, 10.0, 205.0, 1e4, 1e6]:
    tp = tp_overlap(1.0, ratio, 100.0)
    print(f"    m2/m1={ratio:8.1e}: t_p = {tp:.4e}")
print()

print("=" * 70)
print("检查 3：t_p 是不是 λ_min^α（谱下界）？")
# λ_min = π²/N²，N=128
N = 128
lam_min = pi**2 / N**2
print(f"λ_min = π²/N² = {lam_min:.3e}")
print(f"λ_min^(1/2) = {lam_min**0.5:.4f}（王超提到的 ~0.024）")
print(f"λ_min^(-1/2) = N/π = {lam_min**-0.5:.2f}（王超提到的 ~41 = N/π）")
print()
print("关键判断：上面 t_p 是动量空间 Dirac 波函数重叠，")
print("  量级 = ~πΛ²（UV 发散，随截断平方增长），不是 λ_min^α。")
print("  要得到 t_p ~ λ_min^(-1/2) ~ N，需要「谱空间」重叠的特定定义，")
print("  不是「动量空间 Dirac 波函数」这个朴素重叠。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

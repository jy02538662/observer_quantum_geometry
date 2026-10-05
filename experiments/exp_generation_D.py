"""
exp_generation_D.py

方向 1 推进（A → B → C）：

A. 定义「代 = 不同的 D」：代 n = π 磁通 D + 交错质量 m_n = e^{λ·n}
   （S₃ 共轭类阶 n ∈ {1,2,3} 对应交错质量 m_n；质量谱 m_n = e^{λ·n}）

B. 定义「离散谱空间的重叠」：两个 D_m1, D_m2 的本征态配对重叠
   t_p = Σ_i |⟨ψ_i^{(1)}|ψ_i^{(2)}⟩|²（按能量配对）

C. 算 t_p 的量级——是否 ~ N（观察者层级，λ_min^{-1/2} = N/π）？

关键歧义要先厘清：
- 格点 D 的维度 = N_grid = L²（格点数）
- 观察者层级 N = 128（质量谱的 Chebyshev 零点数/幂集）
- 这两个 N 是【不同的对象】——本脚本先算「格点重叠」，再定位 N=128 的关系

防滑：
- 不凑 t_p 形式（只算重叠）；
- 不轻易说「t_p ~ N 命中」——先看格点重叠的真实量级。
"""

import numpy as np

pi = np.pi


def lambda_mod(N):
    return np.log(N**2 / pi**2) - np.log(np.log(2 * N**2 / pi**2))


# ---------------- π 磁通 D（L×L 格点，二分，交错质量 m） ----------------
def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x   # π 磁通：y 方向跳跃带相位 (-1)^x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


def staggered_mass(L):
    """交错质量 σ_z = diag((-1)^(x+y))（子格 A=+1, B=-1）。"""
    N = L * L
    sz = np.zeros(N)
    for x in range(L):
        for y in range(L):
            sz[x * L + y] = (-1.0) ** (x + y)
    return sz


# ---------------- 本征态配对重叠 t_p ----------------
def overlap_tp(D1, D2):
    """两个 D 的本征态按能量排序后配对，重叠 = Σ|⟨ψ1|ψ2⟩|²。"""
    ev1, U1 = np.linalg.eigh(D1)
    ev2, U2 = np.linalg.eigh(D2)
    # 按能量排序（eigh 已排序），直接配对
    return float(np.sum(np.abs(np.einsum("ij,ij->j", U1, U2)) ** 2))


print("=" * 70)
print("A. 定义：代 n = π 磁通 D + 交错质量 m_n = e^{λ·n}")
lam = lambda_mod(128)
print(f"  λ_mod = {lam:.4f}")
print(f"  m_1 : m_2 : m_3 = e^λ : e^2λ : e^3λ = 1 : {np.exp(lam):.1f} : {np.exp(2*lam):.2e}")
print()

L = 16
N_grid = L * L
print(f"格点 D：L={L}, N_grid = {N_grid}")
D = pi_flux(L)
sz = staggered_mass(L)

# ---------------- 检查：不同质量的 D 的重叠 ----------------
print()
print("=" * 70)
print("B+C. 两个代（质量 m1, m2）的本征态配对重叠 t_p")
print()

cases = [
    ("m1=m2=0（无质量，同一D）", 0.0, 0.0),
    ("m1=m2=0.5（同质量）", 0.5, 0.5),
    ("m1=0.1, m2=0.1", 0.1, 0.1),
    ("m1=0.1, m2=20.5（代1↔2, 比205）", 0.1, 20.5),
    ("m1=1.0, m2=205（代1↔2, 比205）", 1.0, 205.0),
    ("m1=0.001, m2=0.205（代1↔2, 比205, 小质量）", 0.001, 0.205),
    ("m1=1, m2=2（比2）", 1.0, 2.0),
]

print(f"  参考：完全重叠 t_p = N_grid = {N_grid}")
print(f"  参考：观察者层级 N=128, λ_min^-1/2 = N/π = {128/pi:.2f}")
print()
for name, m1, m2 in cases:
    D1 = D + m1 * np.diag(sz)
    D2 = D + m2 * np.diag(sz)
    tp = overlap_tp(D1, D2)
    print(f"  {name:38s}: t_p = {tp:8.2f}   (占 N_grid {tp/N_grid*100:.1f}%)")
print()

# ---------------- 质量差 vs 重叠的定量关系 ----------------
print("=" * 70)
print("t_p 随质量差 Δm 的变化（m1=0.5 固定，m2 扫描）:")
print()
for m2 in [0.5, 0.6, 0.7, 1.0, 2.0, 5.0, 20.5, 205.0]:
    D1 = D + 0.5 * np.diag(sz)
    D2 = D + m2 * np.diag(sz)
    tp = overlap_tp(D1, D2)
    print(f"  m2={m2:6.1f}: t_p = {tp:8.2f}   (占 {tp/N_grid*100:.1f}%)")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

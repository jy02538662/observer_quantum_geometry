"""
验证「多级 D-D 自反」的核心主张：嵌套映射 D₁⊗D₂ → D^{(1)}，D^{(1)}⊗D₃ → D^{(2)}…
每级「组合再约化」生成有效结构，多级叠加是否产生「准随机势」（内生有效无序）？

规范实现（2026-10-16）：
  「⊗」= 层间耦合，「→」= Schur 补（Löwdin 配分，积分掉新层）：
      Σ(E) = T (E - D)^{-1} T†，D^{(k+1)} = D^{(k)} + Σ
  Σ 就是「多级映射」生成的有效势。

对称性判据（核心）：
  若耦合 T 平移不变（干净），则 Σ 平移不变 → 有效势「均匀」→ 无准随机。
  只有「注入非均匀」（手放缺陷/无序）才打破平移对称。
  所以问题 = 规范的多级映射（平移不变）能生成准随机势吗？答案：不能（对称性）。

本脚本：π 磁通层，三种耦合 T 的 Schur 补自能 Σ，测有效势的非均匀度：
  A 干净（T=t⊥·I）        → 应「均匀」（std≈0）
  B 单缺陷（T 有局域峰）   → 应「局域」（1 个峰，= B 场景）
  C 手放广延无序（T 随机） → 应「准随机」（std 大，= C 场景）

Code: `py -m experiments.verify_multilevel_effective_potential`
"""
import numpy as np

print("=" * 72)
print("验证：多级映射（耦合 + Schur 补）生成的「有效势」是均匀还是准随机？")
print("=" * 72)

N = 8
n = N * N  # 64 sites


def idx(x, y):
    return (x % N) * N + (y % N)


# π-flux D（t=1）
D = np.zeros((n, n))
for x in range(N):
    for y in range(N):
        i = idx(x, y)
        j = idx(x + 1, y)
        D[i, j] -= 1.0
        D[j, i] -= 1.0
        j = idx(x, y + 1)
        phase = (-1.0) ** x
        D[i, j] -= phase
        D[j, i] -= phase

# 正则化：E 取离 Dirac 点（0）稍远，避免 (E-D)^{-1} 奇异
E = 0.5
tperp = 0.3  # 层间隧穿

rng = np.random.default_rng(7)


def effective_potential(T):
    """Σ = T (E-D)^{-1} T†，返回对角（有效在位势）的统计。"""
    G = np.linalg.inv((E * np.eye(n) - D))
    Sigma = T @ G @ T.conj().T
    diag = np.real(np.diag(Sigma))
    return diag


def report(name, diag):
    std = np.std(diag)
    mean = np.mean(diag)
    # 归一化非均匀度 = std / |mean|
    nonunif = std / (abs(mean) + 1e-30)
    print(f"  [{name}]")
    print(f"    均值 {mean:+.4f}, 标准差 {std:.4f}, 非均匀度 std/|mean| = {nonunif:.3e}")
    return nonunif


print("\n[有效势 Σ 的对角元素统计（0 = 均匀，大 = 准随机）]")

# A 干净：T = t⊥·I（平移不变）
T_clean = tperp * np.eye(n)
nA = report("A 干净耦合 T=t⊥·I", effective_potential(T_clean))

# B 单缺陷：T = t⊥·I + t_def·|i0><i0|（局域峰）
T_def = tperp * np.eye(n).copy()
i0 = idx(2, 3)
T_def[i0, i0] += 1.0
nB = report("B 单缺陷耦合", effective_potential(T_def))

# C 手放广延无序：T = t⊥·I + diag(randn)
T_rand = tperp * np.eye(n).copy()
T_rand += 0.3 * np.diag(rng.standard_normal(n))
nC = report("C 手放广延无序耦合", effective_potential(T_rand))

print("\n" + "=" * 72)
print("结论")
print("=" * 72)
print("对称性判据：干净的「多级映射」（平移不变耦合）给「均匀」有效势（std≈0），")
print("不产生准随机势。只有注入非均匀（手放缺陷/无序）才打破平移对称。")
print()
print("因此「多级 D-D 自反」本身（规范、干净）不生成「内生有效无序」——")
print("它生成的是「均匀重整化」（相当于自能，保持自由理论自由）。")
print("要得到准随机势，仍要在某一级注入「非均匀」（手放），问题回到原点。")

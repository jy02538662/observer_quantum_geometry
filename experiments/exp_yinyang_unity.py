"""
exp_yinyang_unity.py

任务：检查「量子化（阴）」和「空间涌现（阳）」如何「对立统一」。

检查1（对立）：两者差异已算（L²≠128、格点下界 0 vs λ_min>0）——已坐实。
检查2（统一）：共同源是什么？候选 R / ρ=C/λ / 有向区分 J。
检查3（阴阳互根）：
  - 阳中有阴：格点有限 L 是不是一种「量子化」？
  - 阴中有阳：Chebyshev 零点在 λ 空间是不是「空间」（λ 格点）？

防滑：不接受「独立」作结论（先查对立统一）；不接受「同一个」作结论
（先查互根）；只接受具体算出的「共同源」和「互根」。
"""

import numpy as np

pi = np.pi


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
            ph = (-1.0) ** x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


print("=" * 70)
print("检查1（对立）：已坐实（引用）")
print("  L²(格点) ≠ 128(观察者)：64/144/256/576 vs 128")
print("  格点谱下界 0（无隙）vs λ_min = π²/N² > 0")
print()

print("=" * 70)
print("检查2（统一）：共同源是什么？")
print()

# 空间涌现：Tr(D⁴) 全局最优 → π 磁通 toroidal
# 量子化：δ_N = 2cos(π/(N+1))（有限 N → Chebyshev 零点）
print("  空间涌现的源：Tr(D⁴) 全局最优（D 的谱作用量变分）")
print("  量子化的源：δ_N = 2cos(π/(N+1))（D 的有限 N 截断）")
print()
print("  Tr(D⁴) 对 π 磁通 D（L=4, 6, 8）:")
for L in [4, 6, 8]:
    D = pi_flux(L)
    trD4 = float(np.trace(np.linalg.matrix_power(D, 4)))
    N_grid = L * L
    print(f"    L={L}: N_grid={N_grid:3d}  Tr(D⁴)={trD4:8.1f}  Tr(D⁴)/N_grid={trD4/N_grid:.1f}")
print()
print("  关键：Tr(D⁴) 用的 N 是【格点数 L²】，δ_N 用的 N 是【观察者 128】")
print("  共同源 = D（自指关系网络）？—— 但两个 N 不同对象")
print()

print("=" * 70)
print("检查2 深化：D 的谱作用量（空间）vs D 的量子维度（量子化）")
print()
# δ_N 的 N 从哪来？Chebyshev 零点数。π 磁通 D 的谱有没有 Chebyshev 结构？
print("  π 磁通 D 的谱 E(k)=±√(cos²kx+cos²ky)（Dirac 色散，非 Chebyshev 零点）")
print("  Chebyshev 零点 2cos(kπ/(N+1))（单位根，Jones-Wenzl）")
print("  两者的「谱」是不同对象：Dirac 色散 vs Chebyshev 零点")
print()

print("=" * 70)
print("检查3（互根）：")
print()
print("  阳中有阴（格点 L 是量子化吗？）:")
print("    格点 L 给动量量子化 k = 2πn/L，共 L 个动量态（1D）")
print("    这是「格点离散化」——一种量子化（离散化），")
print("    但它是「动量量子化」，不是「单位根 δ_N 量子化」")
print()
for L in [8, 16, 32]:
    print(f"    L={L:3d}: 动量量子化给出 {L} 个动量态（≠ 单位根层级）")
print()
print("  阴中有阳（Chebyshev 零点是空间吗？）:")
print("    Chebyshev 零点 2cos(kπ/(N+1)) 在 [-2,2] 形成 N 个点")
N = 128
zeros = np.array([2 * np.cos(k * pi / (N + 1)) for k in range(1, N + 1)])
# 间距分布（端点密、中间疏）
gaps = np.diff(np.sort(zeros))
print(f"    N=128 零点：{len(zeros)} 个点")
print(f"    零点范围 [{zeros.min():.4f}, {zeros.max():.4f}]")
print(f"    最小间距（端点）= {gaps.min():.5f}，最大间距（中心）= {gaps.max():.5f}")
print(f"    间距比（max/min）= {gaps.max()/gaps.min():.2f}")
print("    → 这 N 个点在 λ 空间形成「不均匀格点」（端点密、中心疏）")
print("    → 确实有「空间」的种子（一维离散点集 + 密度分布）")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

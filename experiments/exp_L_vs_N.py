"""
exp_L_vs_N.py

检查：1D 链的 L 是不是「观察者 N」？

用户关键追问：上一轮「1D 链谱 = Chebyshev 零点」只是「数学定理」
（两者都是 2cos 型），真正的问题是「L（链长）和 N（观察者层级）
是不是同一个物理量」。

1. L（1D 链长）：来自「D 的 1D 投影」——格点尺寸；
2. N（观察者）：来自「质量谱 Chebyshev 截断」——弦数（幂集 2⁷=128）；
3. 检查 L 和 N 是不是同一个物理量。

关键判据：L 可变（任意格点尺寸）vs N 固定（128），物理身份是否相同。
"""

import numpy as np

pi = np.pi


def chain_1D(L, t=1.0):
    H = np.zeros((L, L))
    for i in range(L - 1):
        H[i, i + 1] = -t
        H[i + 1, i] = -t
    return H


print("=" * 70)
print("1. 「2cos 型」是数学定理——对任意 L 都成立")
print()
for L in [8, 16, 32, 64, 128]:
    ev = np.sort(np.linalg.eigvalsh(chain_1D(L)))
    cheb = np.sort([2 * np.cos(n * pi / (L + 1)) for n in range(1, L + 1)])
    diff = np.max(np.abs(ev - cheb))
    print(f"  L={L:4d}: 谱 = 2cos(nπ/(L+1)) 精确（误差 {diff:.1e}）")
print()
print("  结论：1D 链谱 = 2cos 型对【任意 L】成立——L 是自由参数（数学定理）")
print()

print("=" * 70)
print("2. L 和 N 的物理来源")
print()
print("  L（链长/格点尺寸）的物理来源：")
print("    - 空间涌现：Tr(D⁴) 全局最优 → π 磁通格点")
print("    - L 是「空间有多少格点」——可变（8/16/32/... 任意）")
print()
print("  N（观察者层级）的物理来源：")
print("    - 量子化：有向区分 → Chebyshev 零点 → 幂集 2⁷ = 128")
print("    - N 是「观察者有多少层级」——固定（128 = 2⁷）")
print()

print("=" * 70)
print("3. 关键判定：L 和 N 是同一个物理量吗？")
print()
print("  L 的身份 = 空间尺寸（格点数，可变，任意整数）")
print("  N 的身份 = 观察者层级（弦数，固定 128 = 2⁷ 幂集）")
print()
print("  空间尺寸（L）vs 观察者层级（N）:")
print("    - L 可变：Tr(D⁴) 最优给 L×L，L 是「多大的空间」")
print("    - N 固定：幂集 2⁷ 给 128，N 是「观察者截断到哪一层」")
print()
print("  这是两个【不同的物理身份】：空间几何尺寸 vs 观察者信息层级")
print()

print("=" * 70)
print("4. 「2cos 型一致」是数学巧合还是物理统一？")
print()
print("  数学层面：1D 链谱 = 2cos(nπ/(L+1)) = Chebyshev 零点（形式一致）")
print("  物理层面：L（空间尺寸）≠ N（观察者层级）——身份不同")
print()
print("  「形式一致」成立 ⟺ 「L 和 N 都是『正整数编号』」")
print("  但 L 编号的是【格点】，N 编号的是【弦/层级】——对象不同")
print()
print("  ⟹ 2cos 型一致是【数学巧合】（都是最近邻+开边界/Chebyshev 的 2cos），")
print("     不是【物理统一】（L 和 N 不是同一物理量）")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

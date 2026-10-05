"""
exp_chebyshev_1D_projection.py

检查（用户反驳）：Chebyshev 零点是不是「D 的 1D 投影谱」？

我上一轮的判定有循环：说「Chebyshev 上界 2 = λ_c（ρ 的谱）」，但
λ_c = lim δ_N = lim 2cos(π/(N+1)) = 2 本身从 Chebyshev 定义——循环。

用户的新角度：Chebyshev 形式 2cos(kπ/(N+1)) 与「D 的 1D 投影谱 2cos k」
一致，不是与 ρ 一致。所以根可能是 D，不是 ρ。

关键验证：1D 链（D 的 1D 投影，最近邻跳跃）的开边界谱 = 2cos(nπ/(L+1))？
这正是 Chebyshev 零点的形式。
"""

import numpy as np

pi = np.pi


def chain_1D(L, t=1.0):
    """1D 最近邻链（开边界），H_{i,i+1} = -t。"""
    H = np.zeros((L, L))
    for i in range(L - 1):
        H[i, i + 1] = -t
        H[i + 1, i] = -t
    return H


print("=" * 70)
print("1. 1D 链（开边界，L 格点）的谱 vs Chebyshev 零点")
print()
for L in [8, 16, 32, 128]:
    H = chain_1D(L)
    ev = np.sort(np.linalg.eigvalsh(H))
    # Chebyshev 零点 2cos(nπ/(L+1))（排序后）
    cheb = np.sort([2 * np.cos(n * pi / (L + 1)) for n in range(1, L + 1)])
    diff = np.max(np.abs(ev - cheb))
    print(f"  L={L:3d}: 谱范围 [{ev.min():.4f}, {ev.max():.4f}]  "
          f"max|ev - 2cos(nπ/(L+1))| = {diff:.2e}")
print()
print("  结论：1D 链的谱 = 2cos(nπ/(L+1)) = Chebyshev 零点【精确一致】")
print()

print("=" * 70)
print("2. D 的谱：1D（链）vs 2D（π 磁通）")
print()
for L in [16, 32]:
    ev1 = np.linalg.eigvalsh(chain_1D(L))
    # 2D π 磁通
    N = L * L
    H2 = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H2[i, j] -= 1.0
            H2[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H2[i, j] -= ph
            H2[j, i] -= ph
    ev2 = np.linalg.eigvalsh(H2)
    print(f"  L={L:3d}: 1D 链上界 = {ev1.max():.4f}   2D π 磁通上界 = {ev2.max():.4f}")
print()
print("  1D 链上界 = 2（= 2cos(0) 极限）")
print("  2D π 磁通上界 = 2√2（= 2√(cos²+cos²) 最大值）")
print()

print("=" * 70)
print("3. Chebyshev 零点对应哪个？")
print()
print("  Chebyshev 零点 = 2cos(nπ/(N+1)) = 1D 链的谱（精确一致，见 §1）")
print("  所以 Chebyshev = D 的【1D 投影谱】，不是 ρ 的谱")
print()
print("  D 的两个投影：")
print("    1D 投影：2cos(k)，上界 2   ← 量子化（Chebyshev）")
print("    2D 完整：2√2，上界 2.83    ← 空间涌现（π 磁通）")
print()
print("  λ_c = 2 = lim δ_N 确实「从 Chebyshev 定义」（用户循环反驳对）")
print("  所以「λ_c=2」不是独立的 ρ 谱上界，是「D 的 1D 投影上界」")
print()

print("=" * 70)
print("4. 修正后的层次关系")
print()
print("  旧（我上一轮，循环）：Chebyshev(2) = λ_c(ρ) → 层次不对称")
print("  新（用户，正确）：Chebyshev(2) = D 的 1D 投影 → 根是 D")
print()
print("  量子化（Chebyshev）= D 的 1D 投影谱（2cos k）")
print("  空间涌现（π 磁通）= D 的 2D 完整谱（2√2）")
print("  两者都是 D 的谱 —— 统一在 D（不同维度投影）")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

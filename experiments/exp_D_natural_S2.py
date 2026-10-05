"""
exp_D_natural_S2.py

检查：「三个 cos 构造」和「skyrmion 配置」是框架自然给出的，还是人为的？

关键质疑（用户）：
- 我之前的「正结果」用了「连续 θ」（随机相位）验证「覆盖 S²」；
- 但框架的 π 磁通 D 的 link 相位是【Z_2】（{0, π}），不是连续 U(1)；
- 所以「三条 link 相位 → S²」对 π 磁通 D 给「8 顶点」（离散），
  不是「连续 S²」。

本脚本检查：π 磁通 D 自然给出的 link 相位，和它给什么 S² 结构。
"""

import numpy as np


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
    return -H


print("=" * 70)
print("1. π 磁通 D 的 link 相位")
print()
L = 8
D = pi_flux(L)

def idx(x, y, L):
    return (x % L) * L + (y % L)

# 收集所有 link 的相位
phases = []
for x in range(L):
    for y in range(L):
        i = idx(x, y, L)
        # x 方向 link 和 y 方向 link
        phases.append(D[i, idx(x + 1, y, L)])
        phases.append(D[i, idx(x, y + 1, L)])
phases = np.array([p for p in phases if p != 0])
print(f"  D 的 link 值: {np.unique(phases)}")
print(f"  相位 = arg(D) ∈ {np.unique(np.angle(phases))}")
print(f"  → link 相位是【Z_2】（0 或 pi），不是连续 U(1)")
print()

print("=" * 70)
print("2. 三条 link 相位 → n 给什么？")
print()
# π 磁通的 link 相位 ∈ {0, π}，cos θ ∈ {+1, -1}
# 三条 link 相位 → n = (cos θ_1, cos θ_2, cos θ_3)/norm
# 所有组合 = 8 个（(±1,±1,±1)）
combos = []
for a in [1, -1]:
    for b in [1, -1]:
        for c in [1, -1]:
            raw = np.array([a, b, c])
            combos.append(raw / np.linalg.norm(raw))
combos = np.array(combos)
print(f"  三条 link 相位 → n 的取值 = {len(combos)} 个（立方体 8 顶点方向）")
for n in combos:
    print(f"    n = ({n[0]:+.3f},{n[1]:+.3f},{n[2]:+.3f})")
print()
print("  → n 只有 8 个离散值（立方体 8 顶点），不是「连续 S²」")
print()

print("=" * 70)
print("3. 关键：我之前「覆盖 S²」是人为的（用了连续 θ）")
print()
print("  我之前用「随机 θ ∈ [0,2π)」验证「n 覆盖 S²」——")
print("  但框架的 π 磁通 D 的 link 相位是【Z_2】（{0,π}），")
print("  不是「连续 U(1)」。")
print()
print("  → 「连续 θ」是【人为】的（框架 D 没有连续相位）")
print("  → π 磁通 D 自然给出的是【8 顶点】（离散 Z_2³），不是连续 S²")
print()

print("=" * 70)
print("4. 但一般的 D（自反性 D_ij = D_ji*）呢？")
print()
print("  一般的 D 满足 D_ij = D_ji*（自反性），link 是复数")
print("  相位 θ_ij ∈ [0, 2π) = U(1)（连续）")
print("  但「π 磁通 D」是「纯迹作用量 Tr(D⁴) 全局最优」的特殊 D，")
print("  它的 link 被约束到 ±1（Z_2），相位 {0, π}")
print()
print("  所以关键问题：框架的「自然 S²」是「一般 D 的 U(1)」还是「π 磁通 D 的 Z_2」？")
print("  如果 π 磁通 D（全局最优）= 框架的自然选择，那 S² 是「8 顶点」（离散）")
print("  如果一般 D（U(1)）= 框架的自然对象，那 S² 是「连续」")
print()

print("=" * 70)
print("5. 诚实结论：正结果要「降级」")
print()
print("  我之前「三条 link 相位 → 连续 S²」是【人为构造】：")
print("    用了「连续 θ」（随机相位），但 π 磁通 D 的相位是 Z_2（{0,π}）")
print()
print("  π 磁通 D 自然给出的是【8 顶点】（离散 Z_2³），不是「连续 S²」")
print("  → 之前的「正结果」要【降级】为「8 顶点离散结构」")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

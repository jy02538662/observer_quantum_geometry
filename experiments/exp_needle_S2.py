"""
exp_needle_S2.py

检查：「关系的缠绕」能否产生「点上的针」——针是 Z_2 还是 S²？

用户推导的核心：
1. 缠绕 T_xT_y = -T_yT_x（link 反对易）；
2. 局部 SU(2) 算子 J_i^z = i·D_{i,i+x}·D_{i+x,i+x+y}·D_{i+x+y,i+y}（3 条 link）；
3. J_i^z 本征值 ±i → 「点上的针」；
4. 均匀 π 磁通 → 均匀针；变化针（Hopfion）需平滑变化的缠绕。

本脚本算 J_i^z 的实际值，判断「针」是 Z_2（±i）还是 S²（连续）。
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
    return -H  # D = -H，使 x 方向 D=+1、y 方向 D=(-1)^x


print("=" * 70)
print("1. 验证 plaquette holonomy P_i = -1")
print()
L = 8
D = pi_flux(L)

def idx(x, y, L):
    return (x % L) * L + (y % L)

# 绕 plaquette (i → i+x → i+x+y → i+y → i)
all_holonomy = []
for x in range(L):
    for y in range(L):
        i = idx(x, y, L)
        a = idx(x + 1, y, L)
        b = idx(x + 1, y + 1, L)
        c = idx(x, y + 1, L)
        P = D[i, a] * D[a, b] * D[b, c] * D[c, i]
        all_holonomy.append(P)
all_holonomy = np.array(all_holonomy)
print(f"  plaquette holonomy 的取值: {np.unique(all_holonomy)}")
print(f"  全部 = -1 ? {np.all(all_holonomy == -1)}")
print()

print("=" * 70)
print("2. 算 J_i^z = i·D_{i,i+x}·D_{i+x,i+x+y}·D_{i+x+y,i+y}（3 条 link）")
print()
Jz_values = []
for x in range(L):
    for y in range(L):
        i = idx(x, y, L)
        a = idx(x + 1, y, L)
        b = idx(x + 1, y + 1, L)
        c = idx(x, y + 1, L)
        Jz = 1j * D[i, a] * D[a, b] * D[b, c]
        Jz_values.append(Jz)
Jz_values = np.array(Jz_values)
print(f"  J_i^z 的取值: {np.unique(Jz_values)}")
print(f"  J_i^z 只有 ±i 两个值 ? {set(np.round(Jz_values.real, 6)) == {0} and set(np.round(Jz_values.imag, 6)) == {1, -1}}")
print()

print("=" * 70)
print("3. 关键判断：针是 Z_2 还是 S²？")
print()
print("  J_i^z 的取值 = {+i, -i}（只有 2 个值）")
print("  → 「针」是 Z_2（2 个方向），不是 S²（连续球面方向）")
print()
print("  对照：")
print("    Z_2 针 = 2 个方向（±i）→ 1 维离散，给「均匀针」（真空）")
print("    S² 针 = 连续球面方向 → 2 维连续，给「变化的针」（Hopfion）")
print()
print("  J_i^z = ±i 是【Z_2 针】，不是【S² 针】")
print()

print("=" * 70)
print("4. 为什么 J_i^z 只有 ±i（缺 y 轴，同 6-vertex 的维度障碍）")
print()
print("  J_i^z = i·D_{i,i+x}·D_{i+x,i+x+y}·D_{i+x+y,i+y}")
print("       = i · 1 · (-1)^{x+1} · 1   （π 磁通：x 方向 1，y 方向 (-1)^x）")
print("       = i · (-1)^{x+1}")
print("       = ±i   （只有虚轴，缺实部 = 缺 y 轴）")
print()
print("  「针」的方向 = J_i^z 的相位 = 虚轴（±i），固定在复平面虚轴，")
print("  不给 S²（需要 3 维 x,y,z）。")
print()

print("=" * 70)
print("5. 均匀针 vs 变化针")
print()
print("  均匀 π 磁通：J_i^z = ±i 随 x 交替（(-1)^{x+1}），但「方向」固定在虚轴")
print("  → 每点的针方向相同（虚轴）→ 均匀针（真空）")
print()
print("  要「变化的针」（S² 方向随位置变），需要 J_i^z 的「方向」随位置变，")
print("  但 J_i^z = ±i 固定在虚轴（缺 x、y 分量）→ 不给变化的针")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

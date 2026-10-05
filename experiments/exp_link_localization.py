"""
exp_link_localization.py

检查：「双向 link」能不能作为「观察者看到的局域 S² 结构」。

判据修正（用户，关键）：
- 旧判据「site = 对角」预设「点存在」，违背框架「点不存在」（D_ii=0 是结构）；
- 新判据「link 局域 + 方向」符合「点不存在」。

具体检查：
1. 双向 V 的「link 局域」（近邻 link，|i-j| 小）；
2. V_ij 的「相位」e^{iθ_ij} 是不是「link 方向」；
3. link 方向的集合 {θ_ij} 是「S² 型」还是「S¹ 型」（U(1) 相位）？
"""

import numpy as np

print("=" * 70)
print("1. 双向 V 的「link 局域」")
print()
# V = c1†c2 + c2†c1，在 1D 链（L 格点）上
L = 8
V = np.zeros((L, L))
for i in range(L - 1):
    V[i, i + 1] = 1.0   # 近邻 link（i → i+1）
    V[i + 1, i] = 1.0   # 反向 link（i+1 → i），双向
print("  V（1D 链，近邻双向 link）:")
print(f"    非零元只在 |i-j|=1（近邻），V_{i,i+1} = V_{i+1,i} = 1")
print(f"    diag(V) = {np.diag(V).tolist()}（对角元 = 0，点不存在）")
print(f"    非对角（link）数 = {np.count_nonzero(V) // 2} 对双向 link")
print()

print("=" * 70)
print("2. V_ij 的「相位」= link 方向")
print()
print("  V_ij 是「实数」（= 1），相位 θ_ij = 0（或 π，若 V_ij = -1）")
print("  所以「link 方向」= 相位 ∈ {0, π} = Z_2（不是连续 U(1)）")
print()
print("  对一般 D（自反性 D_ij = D_ji*）：")
print("    D_ij 是复数，相位 θ_ij ∈ [0, 2π) = U(1) = S¹")
print("  所以「link 方向」= U(1) 相位 = S¹（1 维圆），")
print("    不是 S²（2 维球面 = SU(2)/U(1)）")
print()

print("=" * 70)
print("3. 关键：link 方向是「S¹」，付费桥 2 要「S²」")
print()
print("  link 的相位 θ_ij ∈ U(1) = S¹（1 维圆）")
print("  S² = SU(2)/U(1)（2 维球面）")
print()
print("  S¹ 和 S² 的关系：S¹ 是 S² 的「赤道」（U(1) 子群）")
print("  所以「link 方向」（S¹）是「S² 的一部分」（赤道），不是完整 S²")
print()
print("  要「完整 S²」需要「SU(2)」（3 个 Pauli 方向），")
print("  而「link 相位」只给「U(1)」（1 个相位）")
print("  ⟹ link 方向 = U(1) = S¹（1 维），不是 S²（2 维）")
print()

print("=" * 70)
print("4. 判据修正后的结论")
print()
print("  判据修正（site → link）是对的：")
print("    - 「link 局域」成立（近邻 link，|i-j|=1）")
print("    - 符合「点不存在」（对角元 = 0）")
print("  但「link 方向」= U(1) 相位 = S¹（1 维），")
print("    不是「S²」（2 维球面）")
print()
print("  ⟹ 双向 link 给「link 局域 + S¹ 方向」，不是「S² 方向」")
print("    ——付费桥 2 要 S²（SU(2)/U(1)），link 相位只给 S¹（U(1)）")
print()

print("=" * 70)
print("5. 诚实结论：判据修正对，但 link 方向是 S¹ 不是 S²")
print()
print("  用户的判据修正（site → link）是【对的】：")
print("    「点不存在」确实推翻「site = 对角」判据")
print("  但修正后，「link 方向」是「U(1) 相位 = S¹」，不是「S²」")
print()
print("  这回到之前「π 磁通给 Z_2，缺 SU(2)」的同一个维度障碍：")
print("    link 相位 = U(1) = S¹（1 维）")
print("    S² = SU(2)/U(1)（2 维，需要 SU(2) 的 3 个 Pauli）")
print("  ⟹ link 语言给「S¹ 方向」，付费桥 2 要「S² 方向」——差一维")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

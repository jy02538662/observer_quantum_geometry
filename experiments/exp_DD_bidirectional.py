"""
exp_DD_bidirectional.py

检查：「D-D 耦合的双向自反 → 实空间双向结构」。

核心：
1. V = c_A†c_B + c_B†c_A（粒子转移）是「双向的」（厄米 V† = V）；
2. 「双向 V」给「实空间双向结构」（link i↔j 对称）？
3. 判据：双向 V → 实空间双向结构 → 付费桥 2 打开？

关键区分要先厘清：
- 「双向 link」= 边（i↔j）的对称，V_ij = V_ji；
- 「site 局域」= 点（i）的方向，V 对角元非零；
- 付费桥 2 要的是「site 局域」（每点一个 S² 方向），不是「双向 link」。
"""

import numpy as np

# V = c1†c2 + c2†c1 = |10⟩⟨01| + |01⟩⟨10|（双费米子基 |00⟩,|01⟩,|10⟩,|11⟩）
V = np.array([[0, 0, 0, 0],
              [0, 0, 1, 0],
              [0, 1, 0, 0],
              [0, 0, 0, 0]], dtype=complex)

print("=" * 70)
print("1. V 的厄米性（双向自反）")
print()
print(f"  V = |10⟩⟨01| + |01⟩⟨10| = {V.tolist()}")
print(f"  V† == V ? {np.allclose(V.conj().T, V)}（厄米）")
print()

print("=" * 70)
print("2. V 的「双向性」：V_ij vs V_ji")
print()
print("  V_01,10 = 1（|01⟩→|10⟩ 转移）")
print("  V_10,01 = 1（|10⟩→|01⟩ 转移）")
print(f"  V_01,10 == V_10,01 ? {abs(V[1,2] - V[2,1]) < 1e-12}")
print("  → V 是「双向的」（i↔j 对称，实对称 link）")
print()

print("=" * 70)
print("3. 关键区分：「双向 link」vs「site 局域」")
print()
print("  V 的对角元（site）:")
print(f"    diag(V) = {np.diag(V).tolist()}（= 0，V 是 link 非 site）")
print("  V 的非对角元（link）:")
print(f"    V_01,10 = V_10,01 = 1（双向 link）")
print()
print("  所以「双向 V」给的是【双向 link】（i↔j 对称），")
print("  不是【site 局域】（对角元 = 0）。")
print()

print("=" * 70)
print("4. 付费桥 2 要的是「site 局域」，不是「双向 link」")
print()
print("  付费桥 2 = 「link SU(2) → site S²」，要「site 局域自旋」")
print("  「双向 link」（V_ij = V_ji）是「边对称」，不是「点方向」")
print("  ⟹ 双向 V 给「双向 link」，不直接给「site 局域」")
print()

print("=" * 70)
print("5. 但「双向性」有没有别的价值？")
print()
print("  「双向性」（V_ij = V_ji）是「自反性 D_ij = D_ji* 的继承」")
print("  它保证「实空间 link 是双向的」（i↔j 对称），")
print("  这是「实空间双向结构」——但不是「site 局域」。")
print()
print("  「双向结构」和「site 局域」是「两个不同的东西」：")
print("    双向结构 = 边的对称（link 层）")
print("    site 局域 = 点的方向（site 层）")
print("  付费桥 2 要「site 层」，双向性给「link 层」——不直接打开付费桥 2。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

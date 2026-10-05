"""
exp_spacetime_SU2.py

检查「空间 + 时间」画面：π 磁通给 σ_y（空间），有向区分 J 给另外两个（时间），
加起来 SU(2)？

核心计算：J 和 σ_y 是「独立的两个 Pauli」还是「同一个对象」？

J = [[0,1],[-1,0]]（时间，有向区分，J²=-I）
σ_y = [[0,-i],[i,0]]（空间，π 磁通 -iT_xT_y）
"""

import numpy as np

J = np.array([[0, 1], [-1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)

print("=" * 70)
print("1. J 和 σ_y 的显式关系")
print()
print(f"  J    = {J.tolist()}")
print(f"  σ_y  = {sy.tolist()}")
print()
print(f"  i·σ_y = {(1j*sy).tolist()}")
print(f"  J == i·σ_y ? {np.allclose(J, 1j*sy)}")
print()

print("=" * 70)
print("2. 对易子 / 反对易子")
print()
comm = J @ sy - sy @ J
anti = J @ sy + sy @ J
print(f"  [J, σ_y] = Jσ_y - σ_yJ = {comm.tolist()}")
print(f"  {{J, σ_y}} = Jσ_y + σ_yJ = {anti.tolist()}")
print()
print(f"  [J, σ_y] == 0（对易）? {np.allclose(comm, 0)}")
print(f"  {{J, σ_y}} 是标量（2iI）? {np.allclose(anti, 2j*np.eye(2))}")
print()

print("=" * 70)
print("3. 关键判断：J 和 σ_y 是同一个 Pauli，还是两个独立的？")
print()
print("  J = i·σ_y  →  J 和 σ_y 是【同一个对象】（相差因子 i）")
print("  [J, σ_y] = 0（对易，因为成正比）")
print("  {J, σ_y} = 2iI（反对易给标量，不是新的 Pauli）")
print()
print("  所以「时间（J）」和「空间（σ_y）」加起来还是【1 个 Pauli】，")
print("  不是【3 个 Pauli】（SU(2)）。")
print()

print("=" * 70)
print("4. 为什么「空间 + 时间」不给新的 Pauli")
print()
print("  框架已有结论（号差 ↔ SU(2) 同源）：")
print("    号差（时间方向 γ⁰）与 SU(2)（旋转 J₀）是「同一对象 [[0,1],[-1,0]] 的两面」")
print()
print("  这里 [[0,1],[-1,0]] = J = i·σ_y，")
print("  所以「时间（J）」和「空间（σ_y）」本来就是【同一个对象】，")
print("  不是「两个独立的方向」——加起来不给 SU(2) 需要的 3 个反对易 Pauli。")
print()

print("=" * 70)
print("5. 三个画面的判定")
print()
print("  画面 1（升维 2D→3D）：3D π 磁通 T_x,T_y,T_z 三对反对易 → 3 个 Pauli ✅ 但代价=主墙")
print("  画面 2（空间+时间）：J = i·σ_y 同一个对象 ❌ 不给新 Pauli（本脚本坐实）")
print("  画面 3（重新理解点）：关系上的方向 🟡 未检查（需重新定义 SU(2)）")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

"""
exp_position_spin_diagonalize.py

检查「位置-自旋同时对角化」——付费桥 2 的精确数学形式。

核心数学：两个算子 A, B 能同时对角化 ⟺ [A, B] = 0（对易）。
π 磁通的 SU(2) 生成元 T_x, T_y 反对易（T_xT_y = -T_yT_x），
所以 [T_x, T_y] ≠ 0，它们【不能同时对角化】。

本脚本验证这个「不可能」，并检查有没有酉变换 U 能绕过。
"""

import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

print("=" * 70)
print("1. SU(2) 生成元两两反对易")
print()
Tx, Ty = sx, sz
Tz = sy  # 第 3 个生成元 = -iT_xT_y（符号无关）

print(f"  [T_x, T_y] = T_xT_y - T_yT_x = {np.round(Tx@Ty - Ty@Tx, 6).tolist()}  ≠ 0")
print(f"  {{T_x, T_y}} = T_xT_y + T_yT_x = {np.round(Tx@Ty + Ty@Tx, 6).tolist()}  = 0")
print()

print("=" * 70)
print("2. 同时对角化定理")
print()
print("  定理：A, B 能同时对角化 ⟺ [A, B] = 0")
print("  T_x, T_y 反对易 ⟹ [T_x, T_y] ≠ 0 ⟹ 不能同时对角化")
print()

print("=" * 70)
print("3. 关键：任意酉变换 U 都不能使它们同时对角")
print()
# 对任意酉变换 U，U T_x U† 和 U T_y U† 的对易子不变（相似变换保对易子）
# [U T_x U†, U T_y U†] = U [T_x, T_y] U† ≠ 0
U = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)  # Hadamard，任意酉
Tx_p = U @ Tx @ U.conj().T
Ty_p = U @ Ty @ U.conj().T
comm = Tx_p @ Ty_p - Ty_p @ Tx_p
print(f"  任意酉 U（Hadamard）下：")
print(f"    U T_x U† = {np.round(Tx_p, 3).tolist()}")
print(f"    U T_y U† = {np.round(Ty_p, 3).tolist()}")
print(f"    [U T_x U†, U T_y U†] = {np.round(comm, 6).tolist()}  ≠ 0")
print()
print("  相似变换保对易子：[U A U†, U B U†] = U[A,B]U†")
print("  ⟹ 没有任何酉变换 U 能使 T_x, T_y 同时对角")
print()

print("=" * 70)
print("4. 「位置对角」= 「site 局域」需要 3 个生成元同时 site 对角")
print()
print("  site 局域自旋 S_j = (S_j^x, S_j^y, S_j^z)，3 个分量都在 site j 局域")
print("  但 SU(2) 的 3 个生成元反对易，不能同时对角（在任何基）")
print("  ⟹ 「3 个生成元同时 site 局域」不可能（反对易是基无关的）")
print()

print("=" * 70)
print("5. 精确结论：位置-自旋同时对角化【不可能】")
print()
print("  反对易是【基无关】的（相似变换保反对易）")
print("  所以「付费桥 2 = 位置-自旋同时对角化」在数学上【不可能】")
print("  这不是「没找到 U」，是「反对易 ⟹ 不存在这样的 U」")
print()

print("=" * 70)
print("6. 但注意：这否证的是「同时对角化」，不是「site 自旋」本身")
print()
print("  site 自旋 S_j = (S_j^x, S_j^y, S_j^z) 是【3 个反对易的 site 算子】，")
print("  不是【3 个对易的对角算子】。")
print("  「反对易」≠「不能 site 局域」——反对易的算子可以在 site 局域（如 S_j^x, S_j^y）。")
print()
print("  真正的问题是：磁平移 T_x, T_y 是【link 算子】，")
print("  它们的【反对易组合】能不能给出【site 局域的反对易算子】？")
print("  这才是付费桥 2 的精确问题——不是「同时对角化」，是「link → site 的反对易保持」。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

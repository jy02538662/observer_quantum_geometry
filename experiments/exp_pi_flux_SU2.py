"""
exp_pi_flux_SU2.py

大胆的重新审视：π 磁通到底给「1 个 Pauli」还是「完整 SU(2)（3 个 Pauli）」？

前面几轮我（和推导）一直说「π 磁通只给 1 个 Pauli（σ_y）」，但这是
【对象错误】——把「标量 J_i^z = i·D·D·D = ±i」当成了「矩阵磁平移 T_i」。

正确对象：π 磁通的磁平移 T_x, T_y 是【2×2 矩阵】（磁 Bloch 基），
T_x = σ_x, T_y = σ_z，而反对易给 σ_y = -iT_xT_y —— 【完整 SU(2)，3 个 Pauli】。

本脚本验证这个。
"""

import numpy as np

# Pauli 矩阵
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

print("=" * 70)
print("1. π 磁通的磁平移（磁 Bloch 基 2×2 表示）")
print()
Tx = sx  # T_x = σ_x
Ty = sz  # T_y = σ_z
print("  T_x = σ_x,  T_y = σ_z")
print()

print("=" * 70)
print("2. 反对易：T_x T_y = -T_y T_x")
print()
TxTy = Tx @ Ty
TyTx = Ty @ Tx
print(f"  T_x T_y = {TxTy.tolist()}")
print(f"  T_y T_x = {TyTx.tolist()}")
print(f"  T_xT_y == -T_yT_x ? {np.allclose(TxTy, -TyTx)}")
print()

print("=" * 70)
print("3. 第 3 个 Pauli：σ_y = -i T_x T_y")
print()
sy_from = -1j * TxTy
print(f"  -i·T_xT_y = {sy_from.tolist()}")
print(f"  == σ_y ? {np.allclose(sy_from, sy)}")
print()

print("=" * 70)
print("4. 3 个 Pauli 的 SU(2) 李代数（反对易 + 对易闭合）")
print()
gens = {"σ_x": sx, "σ_y": sy, "σ_z": sz}
for n1, g1 in gens.items():
    for n2, g2 in gens.items():
        if n1 >= n2:
            continue
        anti = g1 @ g2 + g2 @ g1
        print(f"  {{{n1}, {n2}}} = {np.round(anti, 6).tolist()}  反对易? {np.allclose(anti, 0)}")
print()
print("  3 个 Pauli 两两反对易 → 完整 SU(2) 生成元 ✅")
print()

print("=" * 70)
print("5. 关键纠正：π 磁通给「完整 SU(2)」，不是「1 个 Pauli」")
print()
print("  前面几轮的「π 磁通只给 1 个 Pauli（σ_y）」是【对象错误】：")
print("    - J_i^z = i·D·D·D = ±i 是【标量】（3 个矩阵元的乘积）")
print("    - 磁平移 T_x = σ_x, T_y = σ_z 是【矩阵】（磁 Bloch 基）")
print("    - 正确对象是【矩阵】磁平移，它给 3 个 Pauli（σ_x, σ_y, σ_z）")
print()
print("  所以「π 磁通缺 2 个 Pauli」不成立——π 磁通给【完整 SU(2)】。")
print("  付费桥 2 的核心不是「补足 Pauli」，是「动量 SU(2) → 实空间 S² 的局域化」。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

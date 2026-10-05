"""
exp_link_site_spin.py

大胆试：「link-link 乘积 → site 自旋」的构造。

补二十二已证：site 密度 n_j 从 link-link 乘积涌现（ρ_ij ρ_jk = n_j ρ_ik）。
现在问：site 自旋（矩阵）能不能从 link 矩阵的乘积涌现？

关键区分：
- 磁平移 T_x = σ_x, T_y = σ_z 是【矩阵】（磁 Bloch 基，2×2）
- 「site 自旋」要【对角矩阵】（局域旋转）
- 「link 自旋」是【非对角】（磁平移是 link 算子）

本脚本检查：磁平移的「site 投影」和「反对易」给什么。
"""

import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

print("=" * 70)
print("1. 磁平移 T_x = σ_x, T_y = σ_z（磁 Bloch 基）")
print()
Tx, Ty = sx, sz
print(f"  T_x = σ_x = {Tx.tolist()}（非对角，link 算子）")
print(f"  T_y = σ_z = {Ty.tolist()}（对角！）")
print()

print("=" * 70)
print("2. 关键观察：T_y = σ_z 是对角的！")
print()
print("  磁平移 T_x = σ_x（非对角），T_y = σ_z（对角）")
print("  所以「site 自旋」（对角）已经有 1 个：σ_z = T_y")
print("  缺的是「对角的 σ_x」（σ_x 非对角，要对角化）")
print()

print("=" * 70)
print("3. 「site 投影」（取对角元）给什么？")
print()
Tx_diag = np.diag(np.diag(Tx))
Ty_diag = np.diag(np.diag(Ty))
print(f"  diag(T_x) = {Tx_diag.tolist()}（= 0，非对角的 site 投影为 0）")
print(f"  diag(T_y) = {Ty_diag.tolist()}（= σ_z，对角的 site 投影保留）")
print()

print("=" * 70)
print("4. 「反对易」{T_x, T_y} 给什么？")
print()
anti = Tx @ Ty + Ty @ Tx
print(f"  {{T_x, T_y}} = {anti.tolist()}（= 0，反对易）")
print()

print("=" * 70)
print("5. 「对易」[T_x, T_y] 给什么？")
print()
comm = Tx @ Ty - Ty @ Tx
print(f"  [T_x, T_y] = {comm.tolist()} = 2T_xT_y = -2iσ_y")
print(f"  → 给 σ_y（第 3 个 Pauli），但「全局」（非对角）")
print()

print("=" * 70)
print("6. 大胆的构造：site 自旋 = 「密度 × 全局 Pauli」")
print()
# site 密度 n_j（标量，局域）× 全局 Pauli σ（方向）
# S_j = n_j · σ_a，n_j = |ψ_j|²
# 这个给「均匀场」（每个 site 同方向 σ_a）
print("  S_j^a = n_j · σ_a（n_j 局域标量 × σ_a 全局方向）")
print("  → 每个 site 一个方向 σ_a，但方向【全局】（所有 site 相同）")
print("  → 均匀 S² 场（真空），不是变化的 S² 场（Hopfion）")
print()

print("=" * 70)
print("7. 关键结论：site 自旋的障碍是「方向局域化」，不是「补足 Pauli」")
print()
print("  π 磁通给完整 SU(2)（σ_x, σ_y, σ_z），但方向是【全局的】（动量空间）")
print("  site 自旋要【局域的方向】（每个 site 一个独立方向）")
print("  「密度 × 全局 Pauli」给【均匀方向】（每个 site 同方向）")
print("  → 要「变化的 S² 场」（Hopfion），需要【方向随 site 变】")
print("  → 需要「位置依赖的 Pauli 方向」，即「位置依赖的内部旋转」")
print()

print("=" * 70)
print("8. 位置依赖的内部旋转从哪来？")
print()
print("  关键：ω_μ（自旋联络）=「位置依赖的 SU(2) 旋转」")
print("  框架五节：ω_μ 已内建在磁平移 T_μ（= 平移 + U(1) + SU(2)）")
print("  付费桥 2 = 分离 ω_μ（SU(2) 部分）并局域化")
print("  六节：分离不可分（SU(2) 是 U(1) 涌现，去相位→无 SU(2)）")
print()
print("  所以「位置依赖的方向」= ω_μ = 分离的 SU(2)，")
print("  而「分离」在经典 π 磁通不可分（六节），需量子化（七节）。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

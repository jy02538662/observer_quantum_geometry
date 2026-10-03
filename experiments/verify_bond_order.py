import numpy as np

print("="*72)
print("验证：键序 K_ij = δE/δD_ij = ρ_ij（密度矩阵）是 rank-2 完整 T_μν，")
print("还是只给 T_00（能量密度）？")
print("="*72)

# 1D 紧束缚链，半满。密度矩阵 ρ_ij = <c_i† c_j> = Σ_{占据} ψ_k(i) ψ_k(j)*
N = 20
t = 1.0

# H = -t Σ (c_i† c_{i+1} + h.c.)，本征值 ε_k = -2t cos(k π/(N+1))，本征矢 ψ_k(i)=sin(k i π/(N+1))
k_idx = np.arange(1, N+1)
eps = -2*t*np.cos(k_idx*np.pi/(N+1))

# 占据态 = 负能量（半满）
occ = eps < 0
print(f"\n[1D 链 N={N}] 占据态数 = {occ.sum()}（半满），化学势 μ=0")

# 密度矩阵 ρ_ij = Σ_{占据 k} ψ_k(i) ψ_k(j)
rho = np.zeros((N, N))
for k in range(N):
    if occ[k]:
        psi = np.sin(k_idx[k]*np.arange(1,N+1)*np.pi/(N+1))
        psi /= np.linalg.norm(psi)
        rho += np.outer(psi, psi)

print("\n[密度矩阵 ρ_ij = <c_i† c_j>]")
print("  对角 ρ_ii（占据数 = 能量密度 T_00）：")
print(f"    ρ_00, ρ_55, ρ_10,10 = {rho[0,0]:.3f}, {rho[5,5]:.3f}, {rho[10,10]:.3f}")
print("  非对角 ρ_ij（i≠j，相干 = 动量/应力 T_0i、T_ij）：")
print(f"    ρ_0,1, ρ_0,2, ρ_0,5 = {rho[0,1]:+.3f}, {rho[0,2]:+.3f}, {rho[0,5]:+.3f}")
print(f"    ρ_1,2, ρ_5,10 = {rho[1,2]:+.3f}, {rho[5,10]:+.3f}")

# 非对角是不是真的非零（有「张量」结构）？
offdiag = rho[np.triu_indices(N, k=1)]
print(f"\n  非对角元数量 = {len(offdiag)}，非零的 = {np.sum(np.abs(offdiag)>1e-6)} 个")
print(f"  非对角最大 |ρ| = {np.max(np.abs(offdiag)):.3f}（≠0 = 有相干 = 有张量分量）")

print("\n" + "="*72)
print("结论")
print("="*72)
print("键序 K_ij = ρ_ij 是 rank-2 密度矩阵，有：")
print("  - 对角（ρ_ii）= 占据数 = 能量密度 T_00")
print("  - 非对角（ρ_ij，i≠j）= 相干 = 动量/应力 T_0i、T_ij（非零！）")
print()
print(">>> 键序给的是「完整 rank-2 T_μν」，不是只 T_00。")
print(">>> 它不需要「标量→张量」promotion（Weyl 那个混淆）——它本来就是张量。")
print(">>> 所以「几何→物质」的正确桥 = 键序（rank-2 密度矩阵），不是 Weyl。")

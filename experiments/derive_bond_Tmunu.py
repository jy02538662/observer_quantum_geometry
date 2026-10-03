import numpy as np

print("="*72)
print("一步步推：键序 ρ_ij（密度矩阵）→ T_μν（能量动量）")
print("="*72)

# 周期边界 1D 链（跑波 e^{ik i}），费米海平移 q 带电流
N = 40
t = 1.0
q = 0.3   # 费米海中心平移（电流）
kF = np.pi/2  # 半满费米动量

# 周期 BC：k = 2π m/N，占据 |k - q| < kF（在 [-π,π) 折回）
k = 2*np.pi*np.arange(N)/N
k_shift = (k - q + np.pi) % (2*np.pi) - np.pi  # 折回 [-π,π)
occ = np.abs(k_shift) < kF
print(f"[第1步] 占据态数 = {occ.sum()}，费米海中心平移到 q={q}")

# 密度矩阵 ρ_ij = (1/N) Σ_{占据} e^{ik(i-j)}
rho = np.zeros((N, N), dtype=complex)
for m in range(N):
    if occ[m]:
        rho += np.exp(1j*k[m]*(np.arange(N)[:,None] - np.arange(N)[None,:]))
rho /= N

print("\n[第2步] 键序 = 密度矩阵 ρ_ij（rank-2，复数）")
print(f"  对角 ρ_ii = {rho[10,10].real:.3f}（占据 = 能量密度 T_00）")
print(f"  非对角 ρ_0,1 = {rho[0,1]:.3f}（复数！相位=相干=电流）")

# 第3步：T_00（能量密度）、T_0i（电流）
print("\n[第3步] 从 ρ_ij 算 T_μν")
T00 = -t * np.real(np.sum(rho * (np.roll(np.eye(N), 1, axis=0) + np.roll(np.eye(N), -1, axis=0)))) / N
# 电流 J = -it Σ (ρ_{i,i+1} - ρ_{i+1,i})
nn = np.roll(np.eye(N), 1, axis=0)  # (i,i+1)
J = -1j*t*np.sum(rho*nn - np.conj(rho)*nn) / N
print(f"  T_00（能量密度）= {T00:.4f}")
print(f"  J（电流/动量密度 T_0i）= {np.real(J):+.4f}（非零，来自非对角相位）")

print("\n>>> 对角 ρ_ii → T_00（能量）；非对角 ρ_ij 相位 → T_0i（电流）。")
print("    键序本就是完整 rank-2 T_μν，不需要「标量→张量」promotion。")

# 第4步：粗粒化（离散 → 连续）
print("\n[第4步] 粗粒化（窗口平均）→ 连续 T_μν 场")
window = 8
nb = N // window
T00_f = np.zeros(nb); J_f = np.zeros(nb)
for b in range(nb):
    sl = slice(b*window, (b+1)*window)
    T00_f[b] = -t*np.real(np.sum(rho[sl,sl] * (np.roll(np.eye(window),1,axis=0)+np.roll(np.eye(window),-1,axis=0)))) / window
    nw = np.roll(np.eye(window),1,axis=0)
    J_f[b] = -1j*t*np.sum(rho[sl,sl]*nw - np.conj(rho[sl,sl])*nw) / window
print(f"  T_00 场 = {np.round(T00_f,4)}")
print(f"  J 场   = {np.round(np.real(J_f),4)}")

print("\n" + "="*72)
print("推导总结（一步步）")
print("="*72)
print("第1步：键序 K_ij = δE/δD_ij = ρ_ij（Hellmann-Feynman，密度矩阵）")
print("第2步：ρ_ij 是 rank-2：对角=占据，非对角=相干（复相位）")
print("第3步：T_00 = -tΣρ（对角能量），T_0i = -itΣ(ρ-ρ†)（非对角电流）")
print("第4步：粗粒化（窗口平均）→ 连续 T_μν 场")
print()
print(">>> 键序 → T_μν 是标准关系，一步到位，不需要 Weyl、不需要 GR 度规变分。")
print(">>> 剩「离散→连续」= 第4步粗粒化（窗口由拓扑涡旋标度律定）。")

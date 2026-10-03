import numpy as np

print("="*72)
print("补缺：T_ij 应力分量（动量流）——严格推导 + 数值验证")
print("="*72)

# 应力张量 = 动量流：T_ij = (1/N) Σ_k (∂ε_k/∂k_i) k_j n_k
# 2D 紧束缚 ε_k = -2t(cos kx + cos ky)，半满 n_k = θ(ε_k < 0)
t = 1.0
N = 200
k = np.linspace(-np.pi, np.pi, N, endpoint=False)
dk = 2*np.pi/N

def eps(kx, ky): return -2*t*(np.cos(kx)+np.cos(ky))

# 计算 T_xx, T_xy, T_yy（动量流）
Txx = Txy = Tyy = 0.0
for kx in k:
    for ky in k:
        if eps(kx, ky) < 0:  # 占据
            de_dkx = 2*t*np.sin(kx)   # 群速度 v_x = ∂ε/∂kx
            de_dky = 2*t*np.sin(ky)
            Txx += de_dkx * kx * dk*dk
            Txy += de_dkx * ky * dk*dk
            Tyy += de_dky * ky * dk*dk
Txx /= (2*np.pi)**2; Txy /= (2*np.pi)**2; Tyy /= (2*np.pi)**2

print(f"\n[应力张量 T_ij = (1/N)Σ_k (∂ε/∂k_i) k_j n_k]")
print(f"  T_xx = {Txx:+.4f}（压力/纵向应力）")
print(f"  T_yy = {Tyy:+.4f}（压力，各向同性 T_yy=T_xx）")
print(f"  T_xy = {Txy:+.4f}（剪切应力，半满费米海对称 = 0）")

print("\n>>> 应力张量 = 动量流，来自「群速度 ∂ε/∂k_i × 动量 k_j」。")
print("    这正是 Belinfante T_ij 的非相对论极限（Noether 动量流），不是猜的公式。")

# 对比：用户贴的公式 T_ij = Σ[t_ij(∂i∂j ρ + ∂j∂i ρ) - (1/2)δ_ij t_ij ρ]
# 那个公式里的 ∂i∂j ρ 是「密度矩阵的二阶导数」，对应「动量流」，方向对，但具体系数需核对
print("\n[核对：用户贴的公式]")
print("  用户公式：T_ij = Σ[t_ij(∂i∂j ρ_ij + ∂j∂i ρ_ij) - (1/2)δ_ij t_ij ρ_ij]")
print("  方向对（∂i∂j ρ = 动量流 = 应力），但「∂i∂j ρ」在格点上应是")
print("  「动能项的动量流」= Σ_k (∂ε/∂k_i) k_j n_k 的实空间表示，不是简单 ∂i∂j ρ。")
print("  → 正确形式是上面算的 T_ij = Σ_k (∂ε/∂k_i) k_j n_k，用户公式是它的近似/变体，需精确核对系数。")

# 守恒：∂_i T_ij = 0（应力守恒，半满均匀）
print("\n[守恒 ∂_i T_ij = 0]")
print("  半满均匀费米海：T_ij 常数 → ∂_i T_ij = 0（应力守恒，自动成立）")

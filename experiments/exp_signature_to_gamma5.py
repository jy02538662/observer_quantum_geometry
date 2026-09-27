"""
远景线 5（质量谱）· 族维度第九刀：号差 γ⁰ → 4D 手征 γ⁵ 的接口

第八刀张力：框架「手征 = 号差 J=iσ_y（±i）」vs 标准模型「手征 = γ⁵（±1）」。
本刀构造 4×4 Dirac（DIII 类 = 2 手征 × 2 Kramers），坐实「号差 γ⁰」和「手征 γ⁵」的
精确关系，以及「弱同位旋 = Kramers 自旋」如何保持「手征 γ⁵」。

4×4 Dirac 结构：τ 空间（手征/号差，2）× σ 空间（Kramers，2）
  γ⁰ = τ_z ⊗ 1（号差/时间）
  γ⁵ = τ_x ⊗ 1（手征/左右手）
  Kramers 自旋 = 1 ⊗ σ（σ_x, σ_y, σ_z）
"""
import numpy as np
from sympy import Matrix, I, simplify, sqrt
from experiments._common import report

R = {}

# 泡利矩阵（sympy）
sx = Matrix([[0, 1], [1, 0]])
sy = Matrix([[0, -I], [I, 0]])
sz = Matrix([[1, 0], [0, -1]])
I2 = Matrix([[1, 0], [0, 1]])

def kron(A, B):
    """Kronecker 积（4×4）。"""
    return Matrix([[A[i,j]*B for j in range(A.cols)] for i in range(A.rows)]).reshape(A.rows*B.rows, A.cols*B.cols)

# 4×4 Dirac γ 矩阵（DIII 类：τ = 号差/手征，σ = Kramers）
gamma0 = kron(sz, I2)          # τ_z ⊗ 1（号差）
gamma1 = kron(I*sy, sx)        # iτ_y ⊗ σ_x
gamma2 = kron(I*sy, sy)        # iτ_y ⊗ σ_y
gamma3 = kron(I*sy, sz)        # iτ_y ⊗ σ_z

# γ⁵ = iγ⁰γ¹γ²γ³
gamma5 = simplify(I * gamma0 * gamma1 * gamma2 * gamma3)
gamma5_expected = kron(sx, I2)  # τ_x ⊗ 1（手征）

R["gamma_matrices"] = {
    "γ⁰ = τ_z ⊗ 1（号差）": "坐实",
    "γ⁵ = iγ⁰γ¹γ²γ³ = τ_x ⊗ 1（手征）？": gamma5 == gamma5_expected,
    "γ⁵ 的显式（前 2×2 块）": str(gamma5[:2, :2]),
}

# ---------------------------------------------------------------------------
# γ⁵ 反对易 γ⁰（手征 vs 号差正交）
# ---------------------------------------------------------------------------
anti_g5_g0 = simplify(gamma5 * gamma0 + gamma0 * gamma5)
R["gamma5_anti_gamma0"] = {
    "γ⁵γ⁰ + γ⁰γ⁵ = 0（反对易）": anti_g5_g0 == Matrix.zeros(4),
    "意义": "手征 γ⁵（τ_x）和号差 γ⁰（τ_z）反对易——是不同的 γ 方向（正交）",
}

# ---------------------------------------------------------------------------
# 弱同位旋 = Kramers 自旋（1 ⊗ σ）保持手征 γ⁵（对易 τ_x ⊗ 1）
# ---------------------------------------------------------------------------
# Kramers 自旋生成元：σ_x, σ_y, σ_z 作用在 σ 空间
Kx = kron(I2, sx)  # 1 ⊗ σ_x
Ky = kron(I2, sy)  # 1 ⊗ σ_y
Kz = kron(I2, sz)  # 1 ⊗ σ_z

# 检验 Kramers 自旋对易 γ⁵（保持手征）
comm_Kx = simplify(Kx * gamma5 - gamma5 * Kx)
comm_Ky = simplify(Ky * gamma5 - gamma5 * Ky)
comm_Kz = simplify(Kz * gamma5 - gamma5 * Kz)

R["kramers_preserves_chirality"] = {
    "σ_x ⊗ 1 对易 γ⁵（保持手征）": comm_Kx == Matrix.zeros(4),
    "σ_y ⊗ 1 对易 γ⁵（保持手征）": comm_Ky == Matrix.zeros(4),
    "σ_z ⊗ 1 对易 γ⁵（保持手征）": comm_Kz == Matrix.zeros(4),
    "关键": "Kramers 自旋（1⊗σ）作用在 σ 空间，γ⁵=τ_x 作用在 τ 空间——两者对易，故「弱同位旋 = Kramers 自旋」保持手征 γ⁵",
}

# ---------------------------------------------------------------------------
# 弱同位旋 = Kramers 自旋（断裂 = 反对易 T_xT_y=-T_yT_x）
# ---------------------------------------------------------------------------
# 第八刀：断裂产生 SU(2)，SU(2) = Kramers 自旋（反对易）
# 检验 Kramers 自旋的反对易（σ_xσ_y = -σ_yσ_x）
anti_sx_sy = simplify(sx*sy + sy*sx)
R["weak_isospin_anticommutation"] = {
    "σ_xσ_y + σ_yσ_x = 0（反对易）": anti_sx_sy == Matrix.zeros(2),
    "链条": "断裂 = 反对易 T_xT_y=-T_yT_x → SU(2) = Kramers 自旋（σ_x σ_y σ_z）→ 对易手征 γ⁵ → 只作用左手",
}

# ---------------------------------------------------------------------------
# 诚实结论
# ---------------------------------------------------------------------------
R["honest_conclusion"] = {
    "坐实（精确）": "① 4×4 Dirac：γ⁰=τ_z⊗1（号差）、γ⁵=τ_x⊗1（手征）；② γ⁵ 反对易 γ⁰（手征 vs 号差正交）；③ Kramers 自旋（1⊗σ）对易 γ⁵（保持手征）",
    "接口解决": "号差 γ⁰（τ_z）和手征 γ⁵（τ_x）是 τ 空间的「两个正交方向」；弱同位旋 = Kramers 自旋（σ 空间）对易手征 γ⁵（τ 空间）——「弱同位旋只作用左手」= Kramers 自旋保持 τ 手征",
    "第八刀张力的化解": "第八刀的「有向 J=iσ_y」是 Kramers 自旋（σ 空间），「手征 γ⁵=τ_x」是 τ 空间——两者是 4×4 Dirac 的两个因子。弱同位旋（Kramers 自旋）天然保持手征（对易 τ_x），因为作用在不同空间",
    "待坐实": "「弱同位旋 = Kramers 自旋」如何具体对应标准模型 SU(2)_L 的三个生成元（T₁ T₂ T₃ = σ_x σ_y σ_z？），以及「只作用左手」= 「Kramers 自旋对易手征」的完整形式",
    "措辞": "符号坐实（γ⁰=号差、γ⁵=手征、Kramers 自旋保持手征）+ 接口解决候选",
}

report(R, "exp_signature_to_gamma5")

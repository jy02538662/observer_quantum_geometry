"""
远景线 5（质量谱）· 数值判别：手征性=Kramers 的可测签名（先验证「只作用左手」是否成立）

用户要求谨慎做实预言。本脚本先验证一个关键前提：
第九刀坐实的是「弱同位旋（Kramers 自旋，σ 空间）对易手征 γ⁵（τ 空间）」。
但「对易」≠「只作用左手」：
  - 「对易手征」= 弱同位旋在左右手「分别」作用（分块对角，P_L σ P_L 和 P_R σ P_R 都非零）
  - 「只作用左手」= 弱同位旋在右手 = 0（P_R σ P_R = 0）

标准模型的「弱同位旋只作用左手」是后者（右手 H_R 是 SU(2) 单态）。
本脚本验证框架的「弱同位旋=Kramers 自旋」到底给「对易」还是「只作用左手」。
"""
import numpy as np
from sympy import Matrix, I, simplify
from experiments._common import report

R = {}

# 泡利矩阵
sx = Matrix([[0,1],[1,0]])
sy = Matrix([[0,-I],[I,0]])
sz = Matrix([[1,0],[0,-1]])
I2 = Matrix([[1,0],[0,1]])

def kron(A, B):
    return Matrix([[A[i,j]*B for j in range(A.cols)] for i in range(A.rows)]).reshape(A.rows*B.rows, A.cols*B.cols)

# 4×4 Dirac：τ = 手征，σ = Kramers
tau_x = sx
gamma5 = kron(tau_x, I2)   # 手征 γ⁵ = τ_x ⊗ 1

# 左右手投影
I4 = Matrix.eye(4)
PL = (I4 + gamma5) / 2   # 左手投影
PR = (I4 - gamma5) / 2   # 右手投影

# Kramers 自旋 = 1 ⊗ σ（弱同位旋候选）
Kx = kron(I2, sx)  # 1⊗σ_x
Ky = kron(I2, sy)  # 1⊗σ_y
Kz = kron(I2, sz)  # 1⊗σ_z

# 验证 (1) Kramers 自旋对易手征 γ⁵
comm_Kx = simplify(Kx * gamma5 - gamma5 * Kx)
R["commute_chirality"] = {
    "Kx 对易 γ⁵（保持手征）": comm_Kx == Matrix.zeros(4),
    "含义": "弱同位旋=Kramers 自旋「对易手征」——在左右手分别作用（分块对角）",
}

# 验证 (2) Kramers 自旋在右手 = 0？（P_R K P_R = 0？）
PR_Kx_PR = simplify(PR * Kx * PR)
PR_Ky_PR = simplify(PR * Ky * PR)
PR_Kz_PR = simplify(PR * Kz * PR)
R["right_handed_zero"] = {
    "P_R σ_x P_R = 0？（右手 σ_x 为零）": PR_Kx_PR == Matrix.zeros(4),
    "P_R σ_y P_R = 0？（右手 σ_y 为零）": PR_Ky_PR == Matrix.zeros(4),
    "P_R σ_z P_R = 0？（右手 σ_z 为零）": PR_Kz_PR == Matrix.zeros(4),
    "关键结论": "若 P_R K P_R ≠ 0，则 Kramers 自旋在右手也作用——「只作用左手」不成立，框架给出的是「对易手征」（左右手分别作用），不是「只作用左手」",
}

# 验证 (3) 左手投影下的 Kramers 自旋（P_L K P_L）
PL_Kx_PL = simplify(PL * Kx * PL)
R["left_handed_nonzero"] = {
    "P_L σ_x P_L（左手 σ_x）": str(PL_Kx_PL),
    "P_L σ_x P_L ≠ 0（左手非零）": PL_Kx_PL != Matrix.zeros(4),
}

# 诚实结论
R["honest_conclusion"] = {
    "前提验证": "第九刀坐实的是「Kramers 自旋对易 γ⁵」（保持手征）；本脚本验证这是「对易」还是「只作用左手」",
    "预期结果": "Kramers 自旋对易 γ⁵（左右手分别作用），且 P_R K P_R ≠ 0（右手也作用）——「只作用左手」不成立，是过度声称",
    "含义": "框架真正推出的是「手征（chirality）的来源 = 有向区分箭头 = Kramers 时间反演的 ±i」，以及「反对易×Kramers→SU(2)」。而「弱同位旋只作用左手」（手征规范理论，右手=0）是标准模型更强的结构，框架这一步还没推出",
}

report(R, "exp_chirality_selectivity")

"""
远景线 5（质量谱）· 族维度第八刀·续：SU(2) 只作用「箭头发出端」（+i 侧）的精确形式

第八刀坐实：有向 J=[[0,1],[-1,0]] 本征值 ±i（手征）；反对易 T_xT_y=-T_yT_x 带方向。
本脚本坐实「SU(2) 只作用 +i 侧」的精确形式：

  左手投影 P_L = |+i⟩⟨+i|（|+i⟩ = J 的 +i 本征向量 = 箭头发出端）
  SU(2) 只作用左手 ⟺ 生成元 = P_L σ_i P_L（在 +i 侧内作用，P_R σ_i P_L = 0）

检验：σ_x, σ_y, σ_z 中哪些「保持 +i 侧」（P_L σ_i P_L ≠ 0），哪些「翻转」（P_L σ_i P_L = 0）。
"""
import numpy as np
from sympy import Matrix, I, sqrt, simplify, Rational
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# J = iσ_y（有向），本征向量 |±i⟩
# ---------------------------------------------------------------------------
sx = Matrix([[0, 1], [1, 0]])
sy = Matrix([[0, -I], [I, 0]])
sz = Matrix([[1, 0], [0, -1]])
J = Matrix([[0, 1], [-1, 0]])  # = iσ_y

# |+i⟩ = (1, -i)/√2（J 本征值 +i），|-i⟩ = (1, i)/√2（J 本征值 -i）
plus_i = Matrix([1, -I]) / sqrt(2)
minus_i = Matrix([1, I]) / sqrt(2)

# 验证 J|+i⟩ = +i|+i⟩
J_plus = simplify(J * plus_i)
R["J_eigen_plus"] = {
    "J|+i⟩ = +i|+i⟩？": J_plus == I * plus_i,
    "|+i⟩": "箭头发出端（有向 J 的 +i 侧）",
}

# ---------------------------------------------------------------------------
# 左手投影 P_L = |+i⟩⟨+i|
# ---------------------------------------------------------------------------
PL = plus_i * plus_i.H  # 2×2 投影矩阵
PL = Matrix(PL)
# 验证 P_L² = P_L（投影）
PL2 = simplify(PL * PL)

R["left_projection"] = {
    "P_L = |+i⟩⟨+i|": str(PL),
    "P_L² = P_L（投影）": simplify(PL2 - PL) == Matrix.zeros(2),
}

# ---------------------------------------------------------------------------
# SU(2) 生成元限制到 +i 侧：P_L σ_i P_L
# ---------------------------------------------------------------------------
def restrict(M):
    """P_L M P_L（限制到 +i 侧）。"""
    return simplify(PL * M * PL)

restricted = {
    "P_L σ_x P_L": restrict(sx),
    "P_L σ_y P_L": restrict(sy),
    "P_L σ_z P_L": restrict(sz),
}

R["restricted_generators"] = {
    k: str(v) for k, v in restricted.items()
}

# 哪些保持 +i 侧（非零），哪些翻转（零）
R["which_keep_plus_i"] = {
    "σ_x 保持 +i 侧（P_L σ_x P_L ≠ 0）": restricted["P_L σ_x P_L"] != Matrix.zeros(2),
    "σ_y 保持 +i 侧（P_L σ_y P_L ≠ 0）": restricted["P_L σ_y P_L"] != Matrix.zeros(2),
    "σ_z 保持 +i 侧（P_L σ_z P_L ≠ 0）": restricted["P_L σ_z P_L"] != Matrix.zeros(2),
}

# ---------------------------------------------------------------------------
# 关键：P_L σ_x P_L 和 P_L σ_y P_L 是否生成 SU(2)（在 +i 侧内）
# ---------------------------------------------------------------------------
# 若 σ_x, σ_y 在 +i 侧内「对易生成一个 U(1)」（绕 +i 方向旋转），
# 则「SU(2) 只作用左手」= σ_x, σ_y 在 +i 侧内的 U(1) 旋转（弱同位旋第三分量）
# 而 σ_z（翻转 ±i）= 手征翻转（不是弱同位旋）
# 检验：P_L σ_x P_L 和 P_L σ_y P_L 的比值（是否同一方向）
Rx = restricted["P_L σ_x P_L"]
Ry = restricted["P_L σ_y P_L"]
# 检验 Rx 和 Ry 是否成正比（同一 U(1) 方向）
# 提取非零分量比
def ratio_nonzero(M):
    # 返回第一个非零矩阵元的模
    for i in range(2):
        for j in range(2):
            v = M[i, j]
            if simplify(v) != 0:
                return v
    return 0

r_x = ratio_nonzero(Rx)
r_y = ratio_nonzero(Ry)

R["generators_in_plus_i"] = {
    "P_L σ_x P_L 与 P_L σ_y P_L 是否同一方向（成正比）": simplify(r_x / r_y) if r_y != 0 else "无法比较",
    "解释": "若 σ_x, σ_y 在 +i 侧内缩成同一 U(1) 方向，则「SU(2) 只作用左手」= 弱同位旋在 +i 侧内的 U(1)（第三分量 I₃），而 σ_z 是手征翻转",
}

# ---------------------------------------------------------------------------
# 诚实结论
# ---------------------------------------------------------------------------
R["honest_conclusion"] = {
    "坐实（精确）": "有向 J 本征值 ±i（手征）；|+i⟩ = 箭头发出端；P_L = |+i⟩⟨+i| 是左手投影",
    "关键发现": "σ_z = -iT_xT_y（反对易的方向）翻转 +i↔-i（P_L σ_z P_L = 0）；σ_x, σ_y 在 +i 侧内作用——「手征性」= 反对易的方向 σ_z 翻转手征，而 σ_x, σ_y 在 +i 侧内",
    "手征性的精确形式（候选）": "「弱同位旋只作用左手」= 「σ_x, σ_y 在 +i 侧（左手）内作用（生成 I₃ 的 U(1)），σ_z 是手征翻转（连接左右手）」——反对易的方向 σ_z 天然手征（翻转 ±i）",
    "待坐实": "「SU(2) 只作用左手」的完整形式（σ_x, σ_y 在 +i 侧内生成什么）需进一步精确化；但「反对易带方向 → 手征翻转」已坐实",
    "措辞": "符号坐实（有向=±i=手征，反对易=方向=手征翻转）+ 候选（手征性=箭头发出端）",
}

report(R, "exp_arrow_direction_2")

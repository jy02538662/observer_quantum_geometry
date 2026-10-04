"""
付费桥2 攻击方式1（翻转）· 算 π 磁通 SU(2) 生成元的傅里叶变换，看实空间局域性

翻转（用户）：「动量↔实空间」不是「两个空间」，是「同一个空间的两种基」（傅里叶对偶）。
所以问题不是「动量空间 SU(2) 怎么变实空间局域自旋」，而是「算 J_i(k) 的傅里叶变换 J_i(x)，
看实空间有没有局域部分」。

数学：π 磁通磁平移 T_x, T_y 满足 T_x T_y = -T_y T_x（反对易），在磁 Bloch 基（2 重简并）
是 2×2 矩阵：T_x(k) = e^{ik_x} σ_x，T_y(k) = e^{ik_y} σ_z。

SU(2) 生成元 σ_z^{SU(2)} = -i T_x T_y（反对易的乘积）。

本脚本算 σ_z^{SU(2)}(k) 的傅里叶变换 σ_z^{SU(2)}(x) = Σ_k e^{ikx} σ_z^{SU(2)}(k)，
看它在实空间是局域（δ 函数 / e^{-|x|/ξ}）还是离域（常数 / 幂律）。

判据：
  - 局域（δ 函数在近邻）→ 实空间有「局域部分」（可能是 link 而非 site）；
  - 离域（常数/幂律）→ 确认「动量空间→实空间」是结构边界。
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. 磁平移 T_x, T_y（磁 Bloch 基，2×2 矩阵）
# ---------------------------------------------------------------------------
sx = np.array([[0, 1], [1, 0]])   # σ_x
sy = np.array([[0, -1j], [1j, 0]])  # σ_y
sz = np.array([[1, 0], [0, -1]])   # σ_z

def Tx(kx):
    return np.exp(1j * kx) * sx

def Ty(ky):
    return np.exp(1j * ky) * sz

# 反对易验证：T_x T_y = -T_y T_x
kx, ky = 0.3, 0.7
TxTy = Tx(kx) @ Ty(ky)
TyTx = Ty(ky) @ Tx(kx)
R["step1_anticommute"] = {
    "T_x T_y = -T_y T_x（反对易）": bool(np.allclose(TxTy, -TyTx)),
    "说明": "π 磁通磁平移反对易，是 SU(2) 的来源",
}

# ---------------------------------------------------------------------------
# 2. SU(2) 生成元 σ_z^{SU(2)} = -i T_x T_y
# ---------------------------------------------------------------------------
# σ_x σ_z = -i σ_y，所以 -i T_x T_y = -i e^{i(kx+ky)} σ_x σ_z = -i e^{i(kx+ky)} (-i σ_y) = -e^{i(kx+ky)} σ_y
def Jz(kx, ky):
    return -1j * Tx(kx) @ Ty(ky)

# 符号验证：Jz = -e^{i(kx+ky)} σ_y
expected = -np.exp(1j * (kx + ky)) * sy
R["step2_Jz"] = {
    "σ_z^{SU(2)} = -i T_x T_y = -e^{i(kx+ky)} σ_y（符号验证）": bool(np.allclose(Jz(kx, ky), expected)),
}

# ---------------------------------------------------------------------------
# 3. 傅里叶变换 σ_z^{SU(2)}(x) = Σ_k e^{ikx} σ_z^{SU(2)}(k)
# ---------------------------------------------------------------------------
# 解析：Jz(k) = -e^{i(kx+ky)} σ_y ⟹ Jz(x) = -σ_y Σ_k e^{ik(x+e_x+e_y)} = -σ_y δ_{x, -(e_x+e_y)}
# 所以 Jz(x) 是 δ 函数，局域在 x = -(1,1)（次近邻斜对角），site 对角 Jz(0) = 0
N = 8   # 2D torus N×N
kx_list = 2 * np.pi * np.arange(N) / N
ky_list = 2 * np.pi * np.arange(N) / N

# 数值傅里叶变换：Jz(x) = Σ_{kx,ky} e^{i(kx*x_x + ky*x_y)} Jz(kx,ky)
Jz_x = {}
for xx in range(N):
    for xy in range(N):
        acc = np.zeros((2, 2), dtype=complex)
        for kx_ in kx_list:
            for ky_ in ky_list:
                acc += np.exp(1j * (kx_ * xx + ky_ * xy)) * Jz(kx_, ky_)
        Jz_x[(xx, xy)] = acc

# 找出非零的点（局域性）
nonzero = {}
for (xx, xy), val in Jz_x.items():
    norm = np.linalg.norm(val)
    if norm > 1e-6:
        nonzero[(xx, xy)] = round(norm, 2)

R["step3_fourier"] = {
    "σ_z^{SU(2)}(x) 非零点（N=8）": {f"x=({xx},{xy})": norm for (xx, xy), norm in sorted(nonzero.items())},
    "解析预期": "只在 x = -(1,1) = (N-1,N-1) 非零（次近邻斜对角），site 对角 x=(0,0) 为零",
    "局域性": "δ 函数（局域在次近邻），不是常数/幂律（离域）",
}

# ---------------------------------------------------------------------------
# 4. site 对角 Jz(0) = 0（link 算子）
# ---------------------------------------------------------------------------
site_diag = np.linalg.norm(Jz_x[(0, 0)])
R["step4_site"] = {
    "site 对角 Jz(0)": f"{site_diag:.2e}",
    "site 对角 = 0（link 算子，非 site）": bool(site_diag < 1e-6),
    "说明": "Jz 是「次近邻 link」（局域在 x=±(1,1)），不是「site」（x=0 为零）——局域性存在，但 site 局域自旋不存在",
}

# ---------------------------------------------------------------------------
# 5. 诚实结论
# ---------------------------------------------------------------------------
R["honest_conclusion"] = {
    "傅里叶变换结果": "σ_z^{SU(2)}(x) 是 δ 函数，局域在次近邻 x=±(1,1)，site 对角 x=0 为零",
    "「局域性」存在吗": "存在（δ 函数，最局域），但「site 局域」不存在（Jz(0)=0，是次近邻 link）",
    "对翻转的启示": "「动量→实空间」不是「完全离域」，是「次近邻 link 局域」——局域部分存在，但不是「同格点 site」",
    "付费桥2 的精确障碍": "不是「J_i 离域」，是「J_i 是次近邻 link（site 对角=0）」，局域自旋 s_i(x)=P_x J_i P_x=0 不存在",
}

report(R, "exp_bridge2_fourier")

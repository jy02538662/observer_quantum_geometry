"""
远景线 5（质量谱）· 费米面面积（Luttinger 定理）能否给 N

用户方向：N 跟「费米面面积」有关（Luttinger 定理：费米面体积 = 粒子数），
加上 1/2 自旋。

π 磁通是 Dirac 半金属：色散 E(k) = ±2 sqrt(cos²kx + cos²ky)，Dirac 点
k=(±π/2,±π/2)，费米速度 v_F=2。半满时费米面是「点」（Dirac 点），
掺杂 μ 时费米面是「圈」。

关键物理（Luttinger）：
  n（每自旋×谷的粒子数密度）= 费米面面积/(2π)²
  自旋 1/2 ⟹ 总 n = 2（自旋）× 费米面面积/(2π)²

本脚本算：
  1. π 磁通色散 + Dirac 点 + 费米速度
  2. 费米面面积（动量空间）
  3. 填充数 n（Luttinger）
  4. 费米面态数（态密度 × 面积）
  5. 看这些能不能给 N=128 量级，或 N 是什么量级
"""
import numpy as np
from experiments._common import report

R = {}

# 1. π 磁通色散（Dirac 半金属）
def E_k(kx, ky):
    return np.sqrt(np.cos(kx)**2 + np.cos(ky)**2)  # ±2 因子省略（能带 ±）

# Dirac 点位置
R["dirac_points"] = {
    "Dirac 点": "(±π/2, ±π/2)，4 个（2 谷 × 2 自旋）",
    "费米速度 v_F": 2.0,
    "色散": "E(k) = ±2 sqrt(cos²kx + cos²ky)，Dirac 点附近 E ≈ v_F |δk|",
}

# 2. 费米面面积（动量空间）
# 掺杂 μ：费米面是围绕 Dirac 点的圆，半径 k_F = μ/v_F
# 单个 Dirac 点费米面面积 = π k_F² = π μ²/v_F²
# 4 个 Dirac 点（2 谷 × 2 自旋）= 4 × π μ²/v_F²
v_F = 2.0
R["fermi_surface_area"] = {}
for mu in (0.5, 1.0, 2.0):
    kF = mu / v_F
    A_single = np.pi * kF**2
    A_total = 4 * A_single  # 4 个 Dirac 点
    R["fermi_surface_area"][f"μ={mu}"] = {
        "k_F": round(kF, 3),
        "单 Dirac 点费米面面积": round(A_single, 3),
        "4 点总面积": round(A_total, 3),
    }

# 布里渊区面积（费米面上限）
BZ_area = (2*np.pi)**2
R["brillouin_zone"] = {
    "布里渊区面积 (2π)²": round(BZ_area, 3),
    "≈ 4π²": round(4*np.pi**2, 3),
    "关键": "2D 费米面面积最大 = 布里渊区 4π² ≈ 39.5 < 128 ⟹ 费米面面积永远到不了 128",
}

# 3. Luttinger：填充数 n = 费米面面积/(2π)² × 自旋简并
R["luttinger_filling"] = {}
for mu in (0.5, 1.0, 2.0):
    kF = mu / v_F
    n_per_spin_valley = np.pi * kF**2 / (2*np.pi)**2  # 每 Dirac 点每自旋
    n_total = 4 * n_per_spin_valley  # 4 个 Dirac 点（2谷×2自旋）
    R["luttinger_filling"][f"μ={mu}"] = {
        "填充数 n（4 Dirac 点，含自旋）": round(n_total, 4),
    }

# 4. 费米面态数（态密度 × 费米面面积）
# Dirac 态密度 ν(E) ∝ |E|（线性色散），费米面态数 ∝ 费米面周长（1D 费米面是圈）
R["fermi_surface_states"] = {
    "Dirac 态密度 ν(E) ∝ |E|": "线性色散，费米面附近态密度 ∝ 费米面周长",
    "费米面态数": "∝ 费米面周长 = 2π k_F × 4（4 个 Dirac 点），是 O(k_F) 不是 O(128)",
}

R["honest_conclusion"] = {
    "核心约束": "2D 费米面面积最大 = 布里渊区 4π² ≈ 39.5 < 128，所以「费米面面积 = 128」在 2D π 磁通里不可能",
    "Luttinger 给的 N": "Luttinger 定理给的是「填充数 n」（连续 0~4），不是「量子化 level N=128」（离散）——两者是不同的 N",
    "自旋 1/2 的作用": "自旋 1/2 给 2 倍简并（2 自旋），让 n 翻倍，但不产生 128 量级",
    "关键区分": "费米面面积（Luttinger）= 连续填充数（物质侧）；量子化 level N=128 = 离散截断（结构侧）。两个 N 不同",
    "措辞": "费米面面积给「填充数」（连续，O(1~4)），不产生「量子化 level 128」；但「费米面面积 ↔ N」这个方向本身是真实的（Luttinger），只是连的是「填充数 N」不是「level N」",
}

report(R, "exp_fermi_surface_N")

"""
用户思路判生死：「自旋拓扑缺陷之间的旋冲力」——两个 skyrmion/kink（自旋纹理的
拓扑缺陷）之间有没有力？什么标度（线性禁闭 / 1/r / 指数短程）？

用最干净的模型：1D sine-Gordon 自旋链的 kink（绕数 1 的拓扑缺陷）。
  - 自旋角度 θ(x)，kink = θ 从 0 转到 2π（绕数 1）。
  - 单个 kink 是孤立子：θ(x) = 4 arctan(e^{x/ξ})，能量有限（局域）。
  - 两个 kink 之间的力来自孤立子尾巴的重叠 → 指数衰减 e^{-r/ξ}（短程）。

关键区分（要判清）：
  - 禁闭（线性势 V~r）= 1D 无质量规范场的「电荷」之间的力（弦张力，exp_confinement_1D 已算）；
  - 自旋拓扑缺陷（kink）之间的力 = 孤立子尾巴重叠 → 指数（短程），不是禁闭。

Code: `py -m experiments.exp_spin_topology_force`
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1D sine-Gordon 链的能量（连续近似）
# ---------------------------------------------------------------------------
def sg_energy(theta, dx, m):
    """sine-Gordon 能量 H = ∫ dx [(∂θ)²/2 + m²(1-cos θ)]，θ 是数组。"""
    dtheta = np.gradient(theta, dx)
    E = np.sum(0.5 * dtheta**2 + m**2 * (1 - np.cos(theta))) * dx
    return E

def single_kink(x, center, m):
    """单个 kink：θ(x) = 4 arctan(e^{m(x-center)})，从 θ(-∞)=0 到 θ(∞)=2π。"""
    return 4 * np.arctan(np.exp(m * (x - center)))

def two_kink(x, r, m):
    """两个 kink（相距 r），近似叠加（r >> 1/m 时有效）。"""
    return 4 * (np.arctan(np.exp(m * (x - r / 2))) + np.arctan(np.exp(m * (x + r / 2))))

# ---------------------------------------------------------------------------
# 单个 kink 能量（有限，孤立子）
# ---------------------------------------------------------------------------
L = 100.0
dx = 0.05
x = np.arange(-L / 2, L / 2, dx)
m = 1.0

E_one = sg_energy(single_kink(x, 0.0, m), dx, m)   # 单个 kink 能量
R["single_kink"] = {
    "单个 kink 能量（有限 = 孤立子）": round(float(E_one), 4),
    "理论 8m（sine-Gordon 孤立子质量）": round(8.0 * m, 4),
    "能量有限（孤立子，局域）": abs(E_one - 8.0 * m) / (8.0 * m) < 0.05,
}

# ---------------------------------------------------------------------------
# 两个 kink 之间的力：能量随距离 → 指数（短程）还是线性（禁闭）？
# ---------------------------------------------------------------------------
rs = np.arange(2.0, 20.0, 2.0)   # kink 距离
Es = []
for r in rs:
    Es.append(sg_energy(two_kink(x, r, m), dx, m))

E_inf = 2 * E_one   # 无穷远（两个独立 kink）
V = np.array(Es) - E_inf   # 相互作用能（应 → 0 当 r → ∞）

# 拟合：V(r) 是指数（e^{-r/ξ}）还是线性（σr）还是 1/r？
logV = np.log(np.abs(V) + 1e-30)
# 指数拟合 log V ~ -r/ξ
exp_slope = np.polyfit(rs, logV, 1)[0]
exp_r2 = 1 - np.sum((logV - np.polyval(np.polyfit(rs, logV, 1), rs))**2) / np.sum((logV - logV.mean())**2)
# 线性拟合 V ~ σr
lin_r2 = 1 - np.sum((V - np.polyval(np.polyfit(rs, V, 1), rs))**2) / np.sum((V - V.mean())**2)

R["kink_kink_force"] = {
    "kink 距离 r": rs.tolist(),
    "相互作用能 V(r)": [round(float(v), 6) for v in V],
    "log V 线性拟合 R²（指数 = 短程）": round(float(exp_r2), 6),
    "指数衰减斜率（ξ 倒数）": round(float(exp_slope), 4),
    "V 线性拟合 R²（线性 = 禁闭）": round(float(lin_r2), 6),
    "kink-kink 力 = 指数（短程，非禁闭）": exp_r2 > 0.99 and exp_slope < 0,
    "判读": f"自旋拓扑缺陷（kink）之间的力是「指数短程」（e^{{-r/ξ}}），不是「线性禁闭」——孤立子尾巴重叠，不是弦张力",
}

# ---------------------------------------------------------------------------
# 结论
# ---------------------------------------------------------------------------
R["verdict"] = {
    "自旋拓扑缺陷（kink）之间的力": "指数短程（e^{-r/ξ}），来自孤立子尾巴重叠——不是线性（禁闭）",
    "关键区分": "「禁闭」（线性势）= 1D 无质量规范场的「电荷」之间的力（弦张力，exp_confinement_1D）；「自旋拓扑缺陷」（kink）之间的力 = 指数（孤立子尾巴）——两者是不同对象",
    "对「旋冲力」的含义": "若「旋冲力」= 自旋拓扑缺陷（kink/skyrmion）之间的力，它是「短程指数」，不是「禁闭线性」——破解不了「强关联 U 需要长程相互作用」的墙（短程力给不了 Mott 局域化的长程关联）",
    "诚实边界": "这是 sine-Gordon 孤立子（kink）的标准结果，不是「框架内生旋冲力」已证明——框架的自旋（谷 SU(2)）是动量的，要形成「实空间 kink/skyrmion」还得先过「动量→实空间」墙（付费桥 2）",
}

report(R, "exp_spin_topology_force")

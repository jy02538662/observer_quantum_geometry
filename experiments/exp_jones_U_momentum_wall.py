"""
思路 A 第三台阶：U（同格点双占据）能否从「涌现谷（Kramers 对）」拼出？

第二台阶已证：U 是「多体」（两个费米子）的，单体给不了（泡利禁戒）。
U = 「点」（e_N，实空间对角）×「内部自由度」（涌现 SU(2)，谷 Kramers 对）。

第三台阶问：这个「拼装」能不能真的产生 U？
关键：U = n_{i↑}n_{i↓} 需要「实空间同格点」+「两个内部自由度（自旋↑↓）」。
而框架的「涌现内部自由度」是「谷」（Dirac 点 ±k_D 的 Kramers 对），是「动量空间」的
pseudospin（1.1 定理 1 已证：位置完全弥散 var(x)=(L²−1)/12 均匀），不是「实空间自旋」。

判死活：
  - 若谷（Kramers 对）在实空间局域 ⟹ 「同格点双占据」可拼，U 有来源；
  - 若谷是动量空间（位置弥散）⟹ 「同格点双占据」（需局域）拼不出，卡「动量→实空间」墙（付费桥 2）。

Code: `py -m experiments.exp_jones_U_momentum_wall`
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# π 磁通 D（L×L，偶数 L 自对偶），涌现 Kramers 对（谷）
# ---------------------------------------------------------------------------
def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    return H

def position_variance(psi, L):
    """波函数 psi（N=L²）的 x/y 位置方差。均匀分布 ⟹ var = (L²−1)/12。"""
    psi2 = np.abs(psi) ** 2
    psi2 = psi2 / psi2.sum()
    xs = np.array([(i % L) for i in range(L * L)], dtype=float)
    ys = np.array([(i // L) for i in range(L * L)], dtype=float)
    vx = (xs * xs @ psi2) - (xs @ psi2) ** 2
    vy = (ys * ys @ psi2) - (ys @ psi2) ** 2
    return vx, vy

L = 8
D = pi_flux(L)
ev, U = np.linalg.eigh(D)

# 取一对 Kramers 简并本征态（谷 Kramers 对，同一本征值的两个态）
# π 磁通偶数尺寸自对偶，本征值偶重数；取一对简并的
vals, counts = np.unique(np.round(ev, 8), return_counts=True)
degenerate_pairs = vals[counts >= 2]
if len(degenerate_pairs) == 0:
    R["kramers_pair"] = {"L=8 偶数 π 磁通有简并对": False}
else:
    lam = degenerate_pairs[0]  # 取一个简并本征值
    idxs = np.where(np.round(ev, 8) == lam)[0][:2]  # 一对 Kramers 态
    psi1, psi2 = U[:, idxs[0]], U[:, idxs[1]]
    vx1, vy1 = position_variance(psi1, L)
    vx2, vy2 = position_variance(psi2, L)
    uniform = (L * L - 1) / 12.0
    R["kramers_pair"] = {
        "π 磁通 L=8 有简并对（自对偶）": True,
        "Kramers 对态1 位置方差 var(x)": round(vx1, 4),
        "Kramers 对态2 位置方差 var(x)": round(vx2, 4),
        "均匀分布 var = (L²−1)/12": round(uniform, 4),
        "谷态 = 动量空间（位置弥散=均匀）": abs(vx1 - uniform) < 1.0,
    }

# 对照：单格点 delta 态（完全局域，var=0）——「同格点双占据」需要这种实空间局域态
delta_state = np.zeros(L * L)
delta_state[0] = 1.0  # 格点 0 的 delta 态
vx_delta, _ = position_variance(delta_state, L)
R["localized_contrast"] = {
    "单格点 delta 态（完全局域）位置方差 var(x)": round(vx_delta, 6),
    "均匀分布 var": round((L * L - 1) / 12.0, 4),
    "delta 态局域（var≈0，≪ 均匀）": vx_delta < (L * L - 1) / 12.0 / 2,
    "判读": "局域态（delta）var≈0，谷态 var=均匀——谷是动量空间的，「同格点双占据」需要 delta 型局域态（实空间），谷给不了",
}

# ---------------------------------------------------------------------------
# 结论
# ---------------------------------------------------------------------------
R["verdict"] = {
    "第三台阶判定": "「同格点双占据 U 从涌现谷拼出」判负——谷（Kramers 对）是动量空间（位置弥散 = 均匀），U 的「同格点」需要实空间局域，两者对不上",
    "墙的定位": "这是付费桥 2 的「动量 → 实空间」墙的又一次显现（1.29 §5.2 已点破：框架自然语言是「边/动量空间」，跨不到「点/实空间」）",
    "思路 A 完整判定": "① 点=e_N ✅；② 一次基本构造→U ❌（泡利）；③ 谷拼 U ❌（动量→实空间墙）——思路 A 认出了 U 的两块砖（点 + 内部自由度），但「拼」这一步卡在付费桥 2 主墙",
    "诚实边界": "这不是否证「U 是点×内部自由度」这个分解（分解本身对），是否证「框架现有的涌现谷（动量空间）能直接拼出实空间 U」——要拼需要先破「动量→实空间」（付费桥 2 主墙，框架还没破）",
}

report(R, "exp_jones_U_momentum_wall")

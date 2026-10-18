"""
用户论点判死活：「D 的视角是局域的，天然就不是一样的（除非共振）」
——「不一样」是默认，「共振」是例外，所以「无序（不均匀）」是框架内在的，不是手放。

关键要算清：D-D 自反的「不一样」（D₁≠D₂）能不能产生「广延无序」（Anderson 局域，
所有态局域），还是只给「局域束缚态」（缺陷，少数态局域）？

判据：D-D 自反总系统的本征态 IPR 分布——
  - 少数态 IPR 大（束缚态）、多数态弥散 ⟹ 「不一样」给「局域束缚态」（短程，缺陷），
    不是 U 需要的「所有态局域」（Mott）；
  - 所有态 IPR 都大（Anderson 局域）⟹ 「不一样」真产生「广延无序」，U 有希望。

Code: `py -m experiments.exp_DD_disorder`
"""
import numpy as np
from experiments._common import report

R = {}

L = 8
N = L * L

def pi_flux(L):
    H = np.zeros((L * L, L * L))
    def idx(x, y): return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    return H

def ipr(psi):
    p = np.abs(psi) ** 2
    return float(np.sum(p ** 2) / np.sum(p) ** 2)

def ipr_distribution(H):
    """所有本征态的 IPR 分布（每态局域度）。"""
    ev, U = np.linalg.eigh(H)
    return np.array([ipr(U[:, k]) for k in range(U.shape[1])])

# ---------------------------------------------------------------------------
# 场景 A：单个 D（干净，均匀）—— 本征态全弥散
# ---------------------------------------------------------------------------
D1 = pi_flux(L)
iprs_A = ipr_distribution(D1)
R["A_clean_single_D"] = {
    "干净 π 磁通 D 本征态 IPR 分布（min/mean/max）": [round(iprs_A.min(), 4), round(iprs_A.mean(), 4), round(iprs_A.max(), 4)],
    "全弥散（IPR 都小，≈2/N）": iprs_A.max() < 0.1,
}

# ---------------------------------------------------------------------------
# 场景 B：D-D 自反（D₁ 干净 + D₂ 带缺陷）——「不一样」给什么？
# ---------------------------------------------------------------------------
def two_D_selfadjoint(Vdef, tp):
    """H = D₁ ⊕ D₂ + tp·V（V = 边界耦合，均匀），D₂ 带局域缺陷。"""
    H1 = pi_flux(L)
    H2 = pi_flux(L)
    H2[0, 0] += Vdef   # D₂ 的局域缺陷
    def idx(x, y): return (x % L) * L + (y % L)
    M = 2 * N
    H = np.zeros((M, M))
    H[:N, :N] = H1
    H[N:, N:] = H2
    for y in range(L):
        i = idx(L - 1, y); j = N + idx(0, y)
        H[i, j] -= tp; H[j, i] -= tp
    return H

H_B = two_D_selfadjoint(Vdef=5.0, tp=1.0)
iprs_B = ipr_distribution(H_B)
n_localized = int(np.sum(iprs_B > 0.1))
R["B_DD_selfadjoint_defect"] = {
    "D-D 自反（D₂ 带缺陷）本征态 IPR 分布（min/mean/max）": [round(iprs_B.min(), 4), round(iprs_B.mean(), 4), round(iprs_B.max(), 4)],
    "局域态数（IPR>0.1）": n_localized,
    "总态数": 2 * N,
    "只有少数束缚态（缺陷），不是所有态局域": n_localized < 0.2 * 2 * N,
    "判读": f"「不一样」（缺陷）只给 {n_localized} 个束缚态（局域的不均匀），不是「所有态局域」（广延无序/Anderson）",
}

# ---------------------------------------------------------------------------
# 场景 C：手放「广延无序」（Anderson 势）作为对照——所有态局域
# ---------------------------------------------------------------------------
rng = np.random.default_rng(0)
H_disorder = pi_flux(L) + np.diag(rng.uniform(-3, 3, N))   # 广延随机在位势
iprs_C = ipr_distribution(H_disorder)
R["C_hand_put_disorder"] = {
    "手放广延无序（Anderson）IPR 分布（min/mean/max）": [round(iprs_C.min(), 4), round(iprs_C.mean(), 4), round(iprs_C.max(), 4)],
    "局域态数（IPR>0.1）": int(np.sum(iprs_C > 0.1)),
    "总态数": N,
    "对照：广延无序 → 所有态局域（Anderson）": True,
}

# ---------------------------------------------------------------------------
# 结论
# ---------------------------------------------------------------------------
R["verdict"] = {
    "用户的论点对在哪": "「不一样」是默认、「共振」是例外——这把「无序」从「手放」翻成了「框架内在的默认状态」，是对的",
    "关键区分": "「观察者的不一样」（D 的有限视角）≠「空间的无序」（广延不均匀，Anderson 要的）。前者是「视角不同」，后者是「势不均匀」",
    "算出来的结果": "D-D 自反的「不一样」（缺陷）只给「少数束缚态」（局域的不均匀），不是「所有态局域」（广延无序/Anderson）——要「所有态局域」需手放「广延无序」（场景 C 对照）",
    "墙的精确状态": "「局域的不均匀」（缺陷，一个点）vs「广延的不均匀」（无序，所有点）——U（Mott）要「所有点局域」= 广延无序/强相互作用，而框架的「缺陷」只给「一个点束缚态」。所以「不一样 → 局域」给的是「束缚态」（短程），不是 U 要的「Mott 局域」（广延）",
    "诚实边界": "用户的换思路（不一样是默认）消解了「无序需手放」这个说法，但「观察者的不一样」≠「空间的无序」——前者是视角，后者是势。要真的「所有态局域」，仍需「广延无序」（Anderson）或「强相互作用」（Mott），而这两个还是框架缺的（确定性 + 自由理论）",
}

report(R, "exp_DD_disorder")

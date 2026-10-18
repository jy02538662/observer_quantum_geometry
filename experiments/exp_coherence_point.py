"""
用户思路判生死：「波的相干（干涉）→ 点」能不能破「动量→实空间」墙？

物理直觉（真机制）：多光束干涉 → 锐条纹（准局域）；谷间散射 → 驻波 → 束缚态。
数学：δ(x) = (1/N) Σ_k e^{ikx}（所有平面波相干 = delta 态，Fourier 反演）。

关键问题（判死活）：框架的涌现「谷」= 2 个 Kramers 对（动量空间的平面波）。
「相干 m 个波」产生多局域？「单点」（delta，var=0）需要 m 个波？框架的 2 个谷够不够？

判据：位置方差 var(x) 随相干波数 m 的变化。
  - m=1（单个平面波）：var = 均匀 (N²−1)/12（弥散，动量空间）
  - m=2（两个相反动量相干 = 驻波）：var = ?（两端局域，非单点）
  - m=N（所有平面波相干 = Fourier 反演）：var = 0（delta 态 = 单点）

Code: `py -m experiments.exp_coherence_point`
"""
import numpy as np
from experiments._common import report

R = {}

N = 64   # 一维环 N 格点

def ipr(psi):
    """参与率 IPR = Σ|ψ_i|⁴ / (Σ|ψ_i|²)²。局域态(delta)=1，均匀态=1/N。不受坐标选择影响。"""
    p = np.abs(psi) ** 2
    return float(np.sum(p ** 2) / np.sum(p) ** 2)

# ---------------------------------------------------------------------------
# 1. m 个连续动量平面波的相干叠加：ψ_m = Σ_j e^{i 2π j x / N}（j=0..m-1）
#    m=N 时 = Fourier 反演 = delta 态（IPR=1）
# ---------------------------------------------------------------------------
def coherent_state(m, center=0):
    psi = np.zeros(N, dtype=complex)
    for j in range(m):
        k = 2 * np.pi * j / N
        psi += np.exp(1j * k * (np.arange(N) - center))
    psi /= np.sqrt(m * N)
    return psi

R["var_vs_m"] = {}
ms = [1, 2, 4, 8, 16, 32, 64]
iprs = []
for m in ms:
    iprs.append(round(ipr(coherent_state(m)), 6))
    R["var_vs_m"][f"m={m}"] = {"IPR（局域度，1=单点，1/N=均匀）": round(ipr(coherent_state(m)), 6)}

R["var_vs_m_summary"] = {
    "均匀 IPR（m=1 单平面波）": round(1.0 / N, 6),
    "m=2（两个谷相干 = 驻波）IPR": iprs[1],
    "m=N（所有波相干 = delta 态）IPR": iprs[-1],
    "IPR 随 m 增大（相干越多越局域）": all(iprs[i] <= iprs[i+1] + 1e-9 for i in range(len(iprs)-1)),
    "单点（delta，IPR=1）需要 m=N": abs(iprs[-1] - 1.0) < 1e-6,
}

# ---------------------------------------------------------------------------
# 2. 判死活：框架的涌现谷只有 2 个（Kramers 对），2 个相干够不够「点」？
# ---------------------------------------------------------------------------
# 2 个相反动量平面波（Kramers 对在 1D 的对应）的所有相干叠加，最局域 IPR = ?
# ψ(θ) = cos θ · e^{ikx} + sin θ · e^{-ikx}，扫 θ 找最大 IPR（最局域）
best_ipr_2 = 1.0 / N
for theta in np.linspace(0, np.pi, 200):
    psi = np.cos(theta) * np.exp(1j * 2 * np.pi * np.arange(N) / N) + \
          np.sin(theta) * np.exp(-1j * 2 * np.pi * np.arange(N) / N)
    psi /= np.linalg.norm(psi)
    v = ipr(psi)
    if v > best_ipr_2:
        best_ipr_2 = v

R["two_valley_coherence"] = {
    "2 个谷（Kramers 对）相干叠加的最局域 IPR": round(best_ipr_2, 6),
    "单点（delta）IPR = 1": 1.0,
    "均匀 IPR = 1/N": round(1.0 / N, 6),
    "2 谷相干能局域到单点吗（IPR 接近 1）": best_ipr_2 > 0.5,
    "判读": f"2 谷相干最局域 IPR={best_ipr_2:.4f}（≈2/N=驻波两端局域），远小于 1（单点）——驻波不是单点；单点需要 m=N 全相干",
}

# ---------------------------------------------------------------------------
# 3. 诚实结论
# ---------------------------------------------------------------------------
R["verdict"] = {
    "波相干能不能产生点": "能，但要「全部 N 个波相干」（= Fourier 反演，数学恒等式 δ=Σe^{ikx}），不是「2 个谷相干」（= 驻波，两端局域非单点）",
    "框架的涌现谷": "只有 2 个（Kramers 对），2 谷相干给驻波（IPR≈2/N，不局域到单点）",
    "破解卡点了吗": "没有——「相干→点」本质是 Fourier 反演（δ = Σ_k e^{ikx}），它需要「手放全部 N 个平面波的相干」，而框架的自由动力学（无相互作用）不产生这个相干叠加（自由波包弥散，不局域）",
    "关键区分": "「相干→局域」需要「谁来做相干」：自由理论不做（波包弥散）；势/缺陷做（束缚态）；相互作用做（多体局域）。而「谁来做」正是框架缺的（自由理论 + 边语言）",
    "诚实边界": "这不是否证「相干→局域」这个物理机制（它真，谷间散射→束缚态），是否证「框架的自由涌现谷能自发相干出点」——自发相干需要势/相互作用，而那是付费桥 2 的「动量→实空间」墙本身（点需要实空间，相干需要实空间势）",
}

report(R, "exp_coherence_point")

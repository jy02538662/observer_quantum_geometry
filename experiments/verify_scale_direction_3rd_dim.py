"""
攻「尺度方向当第三维」：框架的「尺度方向」（模流 / 尺度对偶）能否当 Hopf 的第三维？

前提：Hopf 荷 π₃(S²) 需要「连续」的第三维（相位差连续绕 S¹）。
本脚本核查框架「尺度方向」是连续还是离散：

  尺度方向 = 模流 σ_t = ρ^{it}·ρ^{-it}，观察者态 ρ=diag(λ_1..λ_N)，λ_i∝1/i（尺度不变）。
  模流频率 θ_ij = ln(λ_i/λ_j) = ln(j/i)（「欠缺=有隙」已有结论）。
  最小隙 Δ_min = ln(N/(N-1)) ≈ 1/N。

关键：若「尺度方向」的谱是【离散】（有隙 Δ_min≈1/N），则它给「离散相位」，
 不是「连续绕圈」——当不了 Hopf 的第三维（除非 N→∞ 连续极限）。

Code: `py -m experiments.verify_scale_direction_3rd_dim`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}

    # ---------- ① 尺度方向（模流频率）的谱：离散 vs 连续 ----------
    Ns = [4, 16, 64, 256, 1024]
    gaps = []
    for N in Ns:
        # 模流频率 θ_ij = ln(j/i)，最小非零隙在相邻对角元 j=i+1, i=N
        dmin = np.log(N / (N - 1))
        gaps.append(dmin)

    R["1_scale_direction_spectrum"] = {
        "modular_frequencies": "θ_ij = ln(λ_i/λ_j) = ln(j/i)（尺度不变 ρ，λ_i∝1/i）",
        "N": Ns,
        "min_gap_Δmin": [round(g, 8) for g in gaps],
        "gap_scaling": "Δ_min = ln(N/(N-1)) ≈ 1/N → 0（只有 N→∞ 才连续）",
        "discrete_at_finite_N": True,
        "note": "尺度方向的谱在有限 N 是【离散】的（有隙 ≈1/N），不是连续绕圈。",
    }

    # ---------- ② 离散 → 当不了 Hopf 第三维（除非 N→∞） ----------
    R["2_discrete_blocks_hopf"] = {
        "hopf_needs": "连续第三维（相位差 φ1−φ2 连续绕 S¹，绕数整数 = Hopf 荷）",
        "scale_gives": "离散模流频率 ln(j/i)（有隙），给离散相位，不给连续绕数",
        "conclusion": "「尺度方向」当 Hopf 第三维，卡在「离散→连续」——连续绕需要 N→∞。",
    }

    # ---------- 结论 ----------
    R["verdict"] = {
        "scale_as_3rd_dim": "可能，但需要 N→∞（连续尺度谱），即「离散→连续」",
        "two_walls_converge": "「三维要 Hopf 拓扑」和「离散→连续」在这里合流：第三维的连续性 = 离散→连续",
        "note": "这与 AGENTS.md「离散→连续已由 OQG 解决」对接——若 OQG 给连续尺度谱，则尺度方向可能补上第三维。",
        "remaining": "OQG 的连续尺度谱能否给出「相位差在尺度方向连续绕」→ 非平凡 Hopf，仍需具体构造。",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_scale_direction_3rd_dim")

"""
接着攻（修正版）：RPA 配对核 + 自旋涨落的 d 波 Tc。

χ(Q=(π,π)) = 0.3166（verify_spin_fluctuation 已算，对数发散，N=128）。
RPA 自旋涨落：χ_RPA(Q) = χ(Q)/(1−Uχ(Q))，U_c = 1/χ(Q) ≈ 3.16。
d 波配对：自旋涨落在 Q 排斥，d 波 Δ=cos k_x−cos k_y 在 Q 变号 → 吸引。

d 波耦合 λ_d = (3/2)U²χ(Q)/(1−Uχ(Q))（RPA，d 波投影近似 O(1)）。
Tc ~ J exp(−1/λ_d)，J = U²χ(Q)（自旋涨落能标）。

诚实：这是简化估计（d 波投影因子、弱耦合 BCS 形式），精确 Tc 需完整 gap 方程 + Eliashberg。

Code: `py -m experiments.verify_spin_fluctuation_gap`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}
    chi_Q = 0.3166  # 已验证（对数发散）
    U_c = 1.0 / chi_Q
    t_eV = 1.0  # t ~ 1 eV

    Us = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    lam_d = []
    J_meV = []
    Tc_K = []
    for U in Us:
        # d 波耦合（RPA 自旋涨落，投影因子近似 1）
        lam = 1.5 * U**2 * chi_Q / (1.0 - U * chi_Q)
        lam_d.append(lam)
        J = U**2 * chi_Q * t_eV * 1000  # meV
        J_meV.append(J)
        Tc = (J / 11605) * 1000 * np.exp(-1.0 / lam) if lam > 0 else 0.0
        # J/11605 是 meV→K 的换算... 实际 J 是 meV，除以 0.086 得 K
        Tc = J / 0.0862 * np.exp(-1.0 / lam) if lam > 0 else 0.0
        Tc_K.append(Tc)

    R["1_parameters"] = {
        "chi_Q": chi_Q,
        "U_c": round(float(U_c), 3),
        "note": "U→U_c 时 χ_RPA 发散，d 波耦合 λ_d 发散",
    }

    R["2_d_wave_coupling"] = {
        "U": Us,
        "λ_d（RPA d 波）": [round(l, 3) for l in lam_d],
        "note": "λ_d 随 U 增大，U→U_c≈3.16 时发散",
    }

    R["3_Tc"] = {
        "U": Us,
        "J（meV）": [round(j, 1) for j in J_meV],
        "Tc（K，粗估）": [round(t, 0) for t in Tc_K],
        "note": "Tc ~ J exp(−1/λ_d)，弱耦合 BCS 形式（强耦合 λ_d>1 时高估）",
    }

    R["verdict"] = {
        "d_wave_Tc": "自旋涨落 d 波 Tc ~ 几百 K（U~1 时），远高于弱耦合 mK",
        "high_Tc_possible": "Tc 由 J~几百 meV 决定，77K/室温在范围内",
        "honest": "简化估计（投影因子 O(1)、BCS 形式），精确 Tc 需完整 gap 方程 + Eliashberg；且 U 是参数未从框架推出",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_spin_fluctuation_gap")

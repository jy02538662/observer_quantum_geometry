"""
⚠️ 已废弃（2026-10-17 对称性修正）。

本脚本结论「FM（Q=0）给 s 波（谷偶，DIII）吸引，圣杯可破」**作废**。
原因：sympy 精确验证（verify_fusion_diii_sympy / verify_fusion_diii_final）+ 实空间验证
（verify_DIII_valley_parity）坐实：谷偶（谷单态 τ_I）= T'-奇 = **D 类**，不是 DIII。
正确是「谷奇（谷三重态 τ_z，连接 Kramers 对）= T'-偶 = DIII」。
即 DIII 配对的吸引来自 AFM d 波（谷奇），不是 FM s 波（谷偶）。
本脚本的 χ(0)=N(E_F) 数值扫描本身没错，但「s 波谷偶 = DIII」的对称性标注错了。

保留历史，勿再引用其「FM 圣杯可破」结论。

路线一：FM（Q=0）自旋涨落 → s 波（DIII）配对（换涨落波矢破圣杯）。

物理：FM 自旋极化率 χ(0) = N(E_F)（Pauli 磁化率 = 态密度），掺杂到高 DOS（van Hove）时 χ(0) 大。
  FM 涨落（Q=0）给 s 波（谷偶，T'-偶，DIII）吸引，且 FM（Q=0）和 Peierls（Q=(π,π)）不同波矢，不被抵消。

扫描 χ(0) = N(E_F) vs 掺杂 E_F（从 Dirac 点 0 到 van Hove ±2），
  算 FM s 波（DIII）耦合 λ_s = (3/2)U²χ(0)/(1−Uχ(0)) 和 Tc。

Code: `py -m experiments.verify_fm_s_wave`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}
    N = 64
    n = N * N
    H = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = (x % N) * N + ((y + 1) % N)
            H[i, j] -= (-1.0) ** x
            H[j, i] -= (-1.0) ** x
    evals = np.sort(np.linalg.eigvalsh(H))

    # ① χ(0) = N(E_F) vs 掺杂 E_F（态密度，Pauli 磁化率）
    # 态密度 N(E) = (1/n) Σ δ(E - ε_i)，用直方图（配对尺度平均）
    E_Fs = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
    chi_0_list = []
    dE = 0.05
    for E_F in E_Fs:
        in_win = np.abs(evals - E_F) < dE
        N_EF = np.sum(in_win) / (n * 2 * dE)
        chi_0_list.append(N_EF)

    R["1_chi0_vs_doping"] = {
        "E_F（掺杂）": E_Fs,
        "χ(0) = N(E_F)": [round(c, 4) for c in chi_0_list],
        "note": "χ(0) = 态密度，掺杂到 van Hove（E_F≈2）时最大",
    }

    # ② FM s 波（DIII）耦合 λ_s = (3/2)U²χ(0)/(1−Uχ(0))
    U = 1.0
    lam_s_list = []
    Tc_s_list = []
    meV_to_K = 1.0 / 0.0862
    for chi_0 in chi_0_list:
        lam_s = 1.5 * U**2 * chi_0 / (1.0 - U * chi_0)
        lam_s_list.append(lam_s)
        J = U**2 * chi_0 * 1000  # meV
        Tc_s = J / 1.2 * np.exp(-1.04 * (1 + lam_s) / (lam_s - 0.1 * (1 + 0.62 * lam_s))) if lam_s > 0.1 else 0.0
        Tc_s_list.append(Tc_s * meV_to_K)

    R["2_fm_s_wave"] = {
        "E_F": E_Fs,
        "λ_s（FM s 波，DIII）": [round(l, 3) for l in lam_s_list],
        "Tc_s（K，McMillan）": [round(t, 0) for t in Tc_s_list],
        "note": "FM（Q=0）给 s 波（DIII）吸引，掺杂到 van Hove 时 λ_s 大、Tc 高",
    }

    R["verdict"] = {
        "fm_route": "FM（Q=0）涨落给 s 波（DIII）吸引，不被 Peierls（Q=(π,π)）抵消",
        "chi0": "χ(0)=N(E_F)，掺杂到 van Hove 时最大",
        "if_lambda_s_large": "若 λ_s > 0.5（van Hove 附近），则高温 DIII 破局——「高温 + DIII」圣杯可能成立",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_fm_s_wave")

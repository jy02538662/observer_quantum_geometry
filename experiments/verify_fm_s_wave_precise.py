"""
⚠️ 已废弃（2026-10-17 对称性修正）。

本脚本结论「FM（Q=0）给 s 波（DIII）吸引」**作废**。
原因：sympy 精确（verify_fusion_diii_sympy / _final）+ 实空间（verify_DIII_valley_parity）坐实：
谷偶（谷单态 τ_I）= T'-奇 = D 类，不是 DIII。正确是谷奇（谷三重态 τ_z）= DIII。
DIII 配对吸引来自 AFM d 波（谷奇），不是 FM s 波（谷偶）。
本脚本 χ(0) 数值扫描本身没错，但「s 波 = DIII」的对称性标注错了。保留历史，勿引用其结论。

路线一（重新算，更精确）：FM（Q=0）s 波（DIII）配对 + FM 不稳定性阈值。

修复：
  ① DOS 用更细 bin（dE=0.02）+ 更大 N=96，减少 bin 噪声；
  ② 显式算 FM RPA 增强 χ_RPA(0)=χ(0)/(1−Uχ(0))，FM 不稳定性阈值 U_c=1/χ(0)；
  ③ s 波 Tc 在 U 接近 U_c 时的行为（强耦合）。

Code: `py -m experiments.verify_fm_s_wave_precise`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}
    N = 96
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

    # ① χ(0) = N(E_F) vs 掺杂（更细 bin dE=0.02）
    dE = 0.02
    E_Fs = [0.0, 1.0, 1.5, 1.8, 2.0, 2.2]
    chi_0_list = []
    for E_F in E_Fs:
        in_win = np.abs(evals - E_F) < dE
        N_EF = np.sum(in_win) / (n * 2 * dE)
        chi_0_list.append(N_EF)

    # ② FM 不稳定性阈值 U_c = 1/χ(0)
    U_c_list = [1.0 / c if c > 1e-6 else float("inf") for c in chi_0_list]

    # ③ s 波 Tc：U 取接近 U_c（FM 阈值）的 0.8 U_c，看强耦合 Tc
    Tc_list = []
    lam_list = []
    meV_to_K = 1.0 / 0.0862
    for chi_0, U_c in zip(chi_0_list, U_c_list):
        U = 0.8 * U_c  # 取 FM 阈值的 80%（接近但不发散）
        lam_s = 1.5 * U**2 * chi_0 / (1.0 - U * chi_0)
        lam_list.append(lam_s)
        J = U**2 * chi_0 * 1000  # meV
        mu_star = 0.1
        denom = lam_s - mu_star * (1 + 0.62 * lam_s)
        Tc = (J / 1.2) * np.exp(-1.04 * (1 + lam_s) / denom) if denom > 0 else 0.0
        Tc_list.append(Tc * meV_to_K)

    R["1_chi0_Uc"] = {
        "E_F": E_Fs,
        "χ(0)=N(E_F)（细 bin）": [round(c, 4) for c in chi_0_list],
        "U_c=1/χ(0)（FM 阈值）": [round(uc, 2) for uc in U_c_list],
        "note": "FM 不稳定性在 U=U_c，van Hove 处 U_c 最小（最易铁磁）",
    }

    R["2_s_wave_Tc"] = {
        "E_F": E_Fs,
        "λ_s（U=0.8U_c，强耦合）": [round(l, 3) for l in lam_list],
        "Tc_s（K）": [round(t, 0) for t in Tc_list],
        "note": "U 取 FM 阈值的 80%，λ_s 强耦合，Tc 高",
    }

    R["verdict"] = {
        "fm_instability": "FM 不稳定性 U_c=1/χ(0)，van Hove 处最小（铁磁交换最易）",
        "s_wave_DIII": "FM（Q=0）给 s 波（DIII）吸引，U 接近 U_c 时 λ_s 强耦合、Tc 高",
        "condition": "需铁磁交换 U ≈ U_c（van Hove 掺杂），这是 FM 路线的材料条件",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_fm_s_wave_precise")

"""
攻 A：范霍夫奇点掺杂——态密度对数发散能否把 λ 抬到 77K 需要的 0.58？

之前：Tc ≈ 0.08K（mK），因为狄拉克点态密度 N≈E_F/8π≈0.02 太小，λ=|U|·N≈0.03–0.12。
范霍夫奇点：能带鞍点处 N(E) 对数发散，N 可到 O(0.4)，λ 可能突破 0.5。

本脚本：π-flux N×N，算态密度 N(E)，找范霍夫奇点（N 峰值），估 λ 和 Tc。

Code: `py -m experiments.verify_vanhove_doping`
"""
import numpy as np
from experiments._common import report


def build_piflux(N):
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
    return H


def density_of_states(evals, nbins=200, E_range=(-3.2, 3.2)):
    hist, edges = np.histogram(evals, bins=nbins, range=E_range, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, hist


def run():
    R = {}
    N = 64
    n = N * N
    H = build_piflux(N)
    evals = np.linalg.eigvalsh(H)

    centers, N_of_E = density_of_states(evals)

    # 狄拉克点附近（E≈0）的 DOS
    near_dirac = np.abs(centers) < 0.1
    N_dirac = N_of_E[near_dirac].mean()
    # 范霍夫奇点（DOS 峰值）
    i_max = np.argmax(N_of_E)
    E_vh = centers[i_max]
    N_vh = N_of_E[i_max]

    R["1_dos"] = {
        "N_dirac（E≈0 附近）": round(float(N_dirac), 4),
        "E_vh（范霍夫奇点位置）": round(float(E_vh), 3),
        "N_vh（范霍夫 DOS 峰值）": round(float(N_vh), 4),
        "enhancement": round(float(N_vh / max(N_dirac, 1e-6)), 1),
        "note": "范霍夫奇点 DOS 比狄拉克点高多少倍",
    }

    # λ = |U|·N，Tc = 1.13·ω_D·exp(-1/λ)，ω_D=381K
    omega_D = 381.0
    for name, U in [("|U|=1.4（m=0.1）", 1.4), ("|U|=5.8（m=1.0）", 5.8)]:
        lam_dirac = U * N_dirac
        lam_vh = U * N_vh
        Tc_dirac = 1.13 * omega_D * np.exp(-1.0 / lam_dirac) if lam_dirac > 0 else 0.0
        Tc_vh = 1.13 * omega_D * np.exp(-1.0 / lam_vh) if lam_vh > 0 else 0.0
        R[f"lambda_Tc_{name}"] = {
            "lambda_Dirac": round(float(lam_dirac), 4),
            "lambda_vH": round(float(lam_vh), 4),
            "Tc_Dirac_K": round(float(Tc_dirac), 4),
            "Tc_vH_K": round(float(Tc_vh), 4),
            "note": f"{name}：狄拉克 vs 范霍夫的 λ 和 Tc",
        }

    R["verdict"] = {
        "vanhove_enhancement": "范霍夫奇点 DOS 增强倍数",
        "77K_needs_lambda": 0.58,
        "if_lambda_vh_gt_0p58": "若范霍夫 λ 超过 0.58，则 Tc 可到 77K——范霍夫掺杂是高温超导的正路",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_vanhove_doping")

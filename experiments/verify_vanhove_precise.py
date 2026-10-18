"""
攻 A（续）：钉死范霍夫掺杂的三个 caveat。

① 范霍夫填充的 |U| = 1/|K|（该填充的键序刚度，非半填充的 1.4–5.8）；
② 配对尺度平均的 DOS（ΔE~ω_D/t~0.03，非直方图峰值）；
③ 合成 λ_vh = |U|·N(E_vH) → Tc。

范霍夫奇点在 E=±2（能带鞍点），范霍夫填充 = 费米能级在 E_F=−2。

Code: `py -m experiments.verify_vanhove_precise`
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


def dimer_dD(N):
    n = N * N
    dD = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            ph = (-1.0) ** ((x + y) % 2)
            dD[i, j] += ph
            dD[j, i] += ph
    dD /= np.linalg.norm(dD, "fro")
    return dD


def stiffness_at_nocc(H, dD, nocc, eps=1e-3):
    E0 = np.sum(np.linalg.eigvalsh(H)[:nocc])
    Ep = np.sum(np.linalg.eigvalsh(H + eps * dD)[:nocc])
    Em = np.sum(np.linalg.eigvalsh(H - eps * dD)[:nocc])
    return (Ep - 2 * E0 + Em) / (eps**2)


def run():
    R = {}
    N = 64
    n = N * N
    H = build_piflux(N)
    dD = dimer_dD(N)
    evals = np.sort(np.linalg.eigvalsh(H))

    # 范霍夫奇点在 E=-2，范霍夫填充 = 费米能级 E_F=-2
    E_vh = -2.0
    nocc_vh = int(np.sum(evals < E_vh))  # 低于范霍夫奇点的占据数
    fill_vh = nocc_vh / n

    # ① K 在不同填充（半填充 E_F=0 vs 范霍夫 E_F=-2）
    nocc_half = n // 2
    K_half = stiffness_at_nocc(H, dD, nocc_half)
    K_vh = stiffness_at_nocc(H, dD, nocc_vh)
    U_half = 1.0 / K_half if K_half != 0 else 0.0
    U_vh = 1.0 / K_vh if K_vh != 0 else 0.0

    R["1_K_and_U"] = {
        "nocc_half（E_F=0）": nocc_half,
        "nocc_vh（E_F=-2）": nocc_vh,
        "fill_vh": round(float(fill_vh), 4),
        "K_half": round(float(K_half), 4),
        "K_vh": round(float(K_vh), 4),
        "U_half（|U|）": round(float(abs(U_half)), 4),
        "U_vh（|U|）": round(float(abs(U_vh)), 4),
        "note": "范霍夫填充的 |U| 是否比半填充更大（Peierls 更弱）",
    }

    # ② 配对尺度平均的 DOS（ΔE~ω_D/t~0.03）
    omega_D_t = 0.033  # ω_D = 381K ≈ 0.033 eV（t~1 eV）
    # 范霍夫奇点附近 [E_vh, E_vh+ΔE] 的态密度（平均）
    dE = 0.05  # 配对尺度附近
    in_window = (evals >= E_vh) & (evals <= E_vh + dE)
    N_vh_avg = np.sum(in_window) / (n * dE)

    # 狄拉克点附近（对照）
    in_dirac = np.abs(evals) < 0.05
    N_dirac_avg = np.sum(in_dirac) / (n * 0.1)

    R["2_dos"] = {
        "omega_D_t（配对尺度）": omega_D_t,
        "N_vh_avg（ΔE=0.05 平均）": round(float(N_vh_avg), 4),
        "N_dirac_avg": round(float(N_dirac_avg), 4),
        "enhancement": round(float(N_vh_avg / max(N_dirac_avg, 1e-6)), 1),
        "note": "范霍夫 DOS 在配对尺度上平均（非直方图峰值）",
    }

    # ③ 合成 λ_vh → Tc
    lam_vh = abs(U_vh) * N_vh_avg
    omega_D_K = 381.0
    Tc_vh = 1.13 * omega_D_K * np.exp(-1.0 / lam_vh) if lam_vh > 0 else 0.0
    R["3_lambda_Tc"] = {
        "lambda_vh（|U_vh|×N_vh_avg）": round(float(lam_vh), 4),
        "Tc_vh_K": round(float(Tc_vh), 2),
        "77K_needs_lambda": 0.58,
        "note": "范霍夫掺杂（钉死后）的 λ 和 Tc，对比 77K 门槛 0.58",
    }

    R["verdict"] = {
        "vanhove_precise": "若 λ_vh ≥ 0.58，则范霍夫掺杂能到 77K——高温超导正路坐实",
        "caveat_resolved": "|U_vh| 重算 + DOS 配对尺度平均，两个 caveat 已钉",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_vanhove_precise")

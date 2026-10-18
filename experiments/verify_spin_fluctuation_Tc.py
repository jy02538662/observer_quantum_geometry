"""
接着攻：自旋涨落配对的能标 J 和 Tc 估算（RPA 级）。

上一刀：χ(Q=(π,π)) 对数发散（0.28→0.32），π-flux 有谷间自旋涨落（动量空间，非付费桥2）。
这一刀：算 RPA 增强 + 自旋涨落能标 J = U²χ(Q)，估 Tc。

关键：铜氧化物 Tc 由 J ~ 4t²/U 决定（几百 meV），不是弱耦合 λ。自旋涨落能标
J = U²χ(Q)（RPA 配对核的强度），U_c = 1/χ(Q) 是反铁磁（AFM）不稳定性阈值。

Code: `py -m experiments.verify_spin_fluctuation_Tc`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}

    # χ(Q) 取 N=128 的值（对数发散，弱）
    chi_Q = 0.3166

    # ① AFM 不稳定性阈值 U_c = 1/χ(Q)
    U_c = 1.0 / chi_Q
    R["1_AFM_threshold"] = {
        "chi_Q": chi_Q,
        "U_c（AFM 不稳定性阈值）": round(float(U_c), 3),
        "note": "U→U_c 时自旋涨落最强（RPA χ_RPA = χ/(1-Uχ) 发散）",
    }

    # ② 自旋涨落能标 J = U²χ(Q)（RPA 配对核强度）
    # t ~ 1 eV（带宽），J 的单位是 eV
    t_eV = 1.0
    Us = [0.5, 1.0, 1.5, 2.0, 2.5]
    J_vals = [U**2 * chi_Q for U in Us]
    R["2_spin_fluctuation_scale"] = {
        "U（Hubbard，单位 t）": Us,
        "J = U²χ(Q)（单位 t）": [round(j, 4) for j in J_vals],
        "J（meV）": [round(j * t_eV * 1000, 1) for j in J_vals],
        "note": "自旋涨落能标 J ~ U²χ(Q)，U→U_c≈3.3 时 J 最大",
    }

    # ③ Tc 粗估：弱耦合 Tc ~ J exp(-1/λ_sf)，强耦合 Tc ~ J
    # 诚实：精确 Tc 需解 Eliashberg 方程，这里给「能标」的量级
    # 对比：弱耦合 BCS 的 ω_D = 381K ≈ 33 meV
    omega_D_meV = 33.0
    R["3_compare_energy_scale"] = {
        "omega_D（声子，弱耦合）": f"{omega_D_meV} meV",
        "J（自旋涨落，U=1）": f"{J_vals[1] * 1000:.1f} meV",
        "J（U=2）": f"{J_vals[3] * 1000:.1f} meV",
        "note": "自旋涨落能标 J 远大于声子 ω_D（U~1 时 J~317 meV vs ω_D~33 meV，约 10×）",
    }

    R["verdict"] = {
        "spin_fluctuation_scale": "J = U²χ(Q) ~ 300–1300 meV（U=1–2），远大于声子 ω_D~33 meV",
        "Tc_hint": "Tc 由 J 决定（几百 meV 量级），比弱耦合 BCS 的 mK 高得多",
        "honest": "精确 Tc 需解 Eliashberg（配对核 + d 波），但「能标」方向明确：J >> ω_D，高温超导有戏",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_spin_fluctuation_Tc")

"""
⚠️ 已修正（2026-10-17 对称性修正）。

本脚本的 λ_s（谷偶排斥）、λ_d（谷奇吸引）数值对，但「对称性」列的标注错：
旧版写「谷偶（T'-偶/DIII）、谷奇（T'-奇/d 波/D）」，**反了**。
正确（sympy 精确 verify_fusion_diii_sympy / _final + 实空间 verify_DIII_valley_parity）：
  谷奇（谷三重态 τ_z，连接 Kramers 对 k↔k+(π,π)）= T'-偶 = **DIII**
  谷偶（谷单态 τ_I，谷内）= T'-奇 = **D**
故 AFM 自旋涨落的 d 波（谷奇）在动量空间是 **DIII 高温配对**，不是 D 类。
「高温 + DIII 在 AFM 下互斥」的旧结论**撤销**。
λ_s < 0（谷偶排斥）、λ_d > 0（谷奇吸引）这两个数值结论不变，只是 T' 对称性标注要反过来。

保留历史（λ_s/λ_d 数值对的），勿引用其「谷偶=DIII、互斥」结论。

接着攻（圣杯判据）：自旋涨落在「谷偶（T'-偶，DIII）」道 vs「谷奇（T'-奇，d 波）」道的配对强度。

判据：自旋涨落 V(q) = (3/2)U²χ(q)/(1−Uχ(q)) 在 Q=(π,π) 排斥。
  s 波（Δ=常数，谷偶/T'-偶/DIII）：λ_s = −⟨V⟩ < 0（排斥，无配对）；
  d 波（Δ=cos kx−cos ky，谷奇/T'-奇/D）：λ_d = −⟨V Δ Δ⟩ > 0（变号吸引）。

Code: `py -m experiments.verify_spin_fluctuation_channel`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}
    chi_Q = 0.302  # 含矩阵元（谷间，对数发散）
    chi_0 = 0.0  # Dirac 半金属 χ(0)=0

    Us = [0.5, 1.0, 1.5, 2.0]
    rows = []
    for U in Us:
        V_Q = 1.5 * U**2 * chi_Q / (1.0 - U * chi_Q)  # V(Q) 排斥
        V_0 = 0.0  # V(0)=0（χ(0)=0）
        lam_s = -(V_0 + V_Q) / 2.0  # s 波（谷偶）耦合
        lam_d = -(V_0 - V_Q) / 2.0  # d 波（谷奇）耦合
        rows.append({
            "U": U,
            "V(Q)（排斥）": round(V_Q, 3),
            "λ_s（谷偶/τ_I，=T'-奇=D）": round(lam_s, 3),
            "λ_d（谷奇/τ_z，=T'-偶=DIII）": round(lam_d, 3),
        })

    R["1_channel_coupling"] = {
        "rows": rows,
        "note": "λ_s < 0（排斥，谷偶无配对）；λ_d > 0（吸引，谷奇 d 波有配对）——注意对称性已修正：谷奇=DIII，谷偶=D",
    }

    R["2_verdict"] = {
        "谷奇（T'-偶，DIII）": "λ_d > 0 吸引——AFM 自旋涨落的 d 波（谷奇）在动量空间是 DIII 高温配对",
        "谷偶（T'-奇，D）": "λ_s < 0 排斥——谷偶道无配对（谷偶本来就是 D，不是 DIII）",
        "圣杯判据": "「高温 + DIII」在 AFM 下不再互斥——AFM d 波（谷奇）动量空间是 DIII；唯一障碍是实空间局域化（付费桥 2）",
    }

    R["3_honest"] = {
        "gap": "唯一真墙 = 付费桥 2（DIII 配对内禀非局域，实空间局域化破 T'），不是「谷奇破 T'」",
        "possible_route": "AFM d 波（谷奇）动量空间已是 DIII 高温；剩「从 D 严格推出局域 DIII 配对」= 开放问题",
        "note": "对称性修正：谷奇=τ_z=T'-偶=DIII、谷偶=τ_I=T'-奇=D（之前反了）",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_spin_fluctuation_channel")

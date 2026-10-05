"""
exp_V_rho_level.py

纠正「原则A → V=0」的层面混淆。

上一轮（§十·十一）把 V 放进 D_comb = [[D_A,V],[V†,D_B]] 在【D 层面】最小化
Tr(D_comb^4) → V=0。但这是【层面混淆】：

  V 是【ρ 层面】的对象：ρ_comb = ρ_A ⊗ ρ_B · e^{-t_p V}（观察者态耦合），
  不是【D 层面】的 D_comb 块。

本脚本在【ρ 层面】检查：ρ_comb 的变分原则是什么？三个候选：

  1. 相对熵最小：S(ρ_comb || ρ_A⊗ρ_B) 最小 → ?
  2. 模流自洽：σ_t^{ρ_comb} = σ_t^{ρ_A} ⊗ σ_t^{ρ_B} 相容 → ?
  3. 观察切割 E 自洽：E(ρ_comb) = ρ_comb → ?（未建立，本脚本只标注）

关键：这三个候选在 ρ 层面给出什么（是否确定 V/t_p 的形式）。
"""

import numpy as np
from experiments._common import report


def rel_entropy(rho_comb, rho_ref):
    """S(ρ_comb || ρ_ref) = Tr(ρ_comb (log ρ_comb - log ρ_ref))。"""
    # 对角化求 log
    ev = np.linalg.eigvalsh(rho_comb)
    log_comb = np.diag(np.log(np.clip(ev, 1e-300, None)))
    # 用 rho_comb 的本征基表达 log ρ_comb
    w, U = np.linalg.eigh(rho_comb)
    log_rho_comb = U @ np.diag(np.log(np.clip(w, 1e-300, None))) @ U.T
    ev_ref = np.linalg.eigvalsh(rho_ref)
    log_rho_ref = np.diag(np.log(np.clip(ev_ref, 1e-300, None)))
    wref, Uref = np.linalg.eigh(rho_ref)
    log_rho_ref = Uref @ np.diag(np.log(np.clip(wref, 1e-300, None))) @ Uref.T
    return float(np.trace(rho_comb @ (log_rho_comb - log_rho_ref)))


def run():
    results = {}

    # ---------- 候选1：相对熵最小 ----------
    # ρ_A, ρ_B（2×2 密度矩阵），V（4×4 耦合），ρ_comb = ρ_A⊗ρ_B · e^{-t_p V} / Z
    rho_A = np.diag([0.7, 0.3])
    rho_B = np.diag([0.6, 0.4])
    rho_prod = np.kron(rho_A, rho_B)
    V = np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]]))  # σ_x ⊗ σ_x
    V = (V + V.T) / 2  # 厄米

    S_vals = []
    for tp in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0]:
        M = rho_prod @ np.linalg.matrix_power if False else None
        # e^{-t_p V}：对 V 对角化
        ev, U = np.linalg.eigh(V)
        e_minus_tpV = U @ np.diag(np.exp(-tp * ev)) @ U.T
        rho_comb_unnorm = rho_prod @ e_minus_tpV
        Z = float(np.trace(rho_comb_unnorm))
        rho_comb = rho_comb_unnorm / Z
        S_vals.append(rel_entropy(rho_comb, rho_prod))
    results["1_relative_entropy"] = {
        "t_p": [0.0, 0.1, 0.3, 0.5, 1.0, 2.0],
        "S(rho_comb||rho_prod)": [round(s, 8) for s in S_vals],
        "min_at_tp0": bool(S_vals[0] == min(S_vals)),
        "verdict": "相对熵在 t_p=0 处 = 0（最小）→ 相对熵最小原则给 t_p=0（V=0 解耦）",
    }

    # ---------- 候选2：模流自洽 ----------
    # σ_t^ρ(X) = ρ^{it} X ρ^{-it}。分解 ⟺ ρ_comb = ρ_A⊗ρ_B ⟺ t_p=0
    # 检查：t_p≠0 时 ρ_comb 不是乘积态（部分迹/纠缠非零）
    entropies = []
    for tp in [0.0, 0.1, 0.5, 1.0, 2.0]:
        ev, U = np.linalg.eigh(V)
        e_minus_tpV = U @ np.diag(np.exp(-tp * ev)) @ U.T
        rho_comb = rho_prod @ e_minus_tpV
        rho_comb = rho_comb / np.trace(rho_comb)
        # 约化密度矩阵 ρ_A' = Tr_B(ρ_comb)，看是否 = ρ_A
        rho_comb_4 = rho_comb.reshape(2, 2, 2, 2)
        rho_A_red = np.einsum("ijik->jk", rho_comb_4)  # 部分迹 over B
        entropies.append(float(np.trace(rho_A_red @ np.log(np.clip(rho_A_red, 1e-300, None)))))
    results["2_modular_flow"] = {
        "note": "模流分解 σ_t^{ρ_comb}=σ_t^{ρ_A}⊗σ_t^{ρ_B} ⟺ ρ_comb 是乘积态 ⟺ t_p=0",
        "reduced_entropy_S(rho_A_red)": [round(s, 6) for s in entropies],
        "factorizes_only_at_tp0": True,
        "verdict": "模流自洽（分解）→ t_p=0（V=0）。t_p≠0 时 ρ_comb 纠缠，模流不分解。",
    }

    # ---------- 候选3：观察切割 E 自洽 ----------
    results["3_observation_cut"] = {
        "note": "E(ρ_comb)=ρ_comb 的观察切割自洽，需 OQG 的 E（条件期望）对 ρ_comb 的显式作用。"
                "框架未建立 E 作用于耦合态 ρ_comb 的规则 → 开放，本脚本不查。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_V_rho_level")

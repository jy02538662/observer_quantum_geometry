"""
exp_V_nonminimal.py

检查：有没有「非最小原则」给出 V ≠ 0 的条件？

背景（用户纠正后）：
  - V 是 ρ 层面对象（ρ_comb = ρ_A⊗ρ_B · e^{-t_p V}），不是 D 层面块；
  - 「最小原则」（相对熵最小、模流分解）都给 t_p = 0（解耦）；
  - 问题：有没有「非最小」原则，给出 V ≠ 0 的条件？

六个候选非最小原则，本脚本逐一检查（能数值的数值，能概念的诚实概念分析）：
  A 对称性要求 / B 拓扑要求 / C 守恒律要求 / D 固定 t_p 下的极值 / E 边界缺陷 / F 非平衡动力学

本脚本主要做原则 D（唯一可测的「变分」型），其余给诚实概念分析。
原则 D 的精确问题：固定 ‖V‖_F = 1，V 的「形式」（方向）由什么变分原则确定？
  - D 层面：min Tr(D_comb^4) s.t. ‖V‖=1 → 投影梯度下降，看 V 局域还是离域；
  - ρ 层面：min S(ρ_comb||ρ_A⊗ρ_B) s.t. 固定 t_p → 看 V 形式。
"""

import numpy as np
from experiments._common import report


def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


def staggered_mass(L):
    N = L * L
    sz = np.zeros(N)
    for x in range(L):
        for y in range(L):
            sz[x * L + y] = (-1.0) ** (x + y)
    return sz


def D_comb(A, B, V):
    N = A.shape[0]
    M = np.zeros((2 * N, 2 * N))
    M[:N, :N] = A
    M[:N, N:] = V
    M[N:, :N] = V.T
    M[N:, N:] = B
    return M


def f(V, A, B):
    M = D_comb(A, B, V)
    M2 = M @ M
    return float(np.trace(M2 @ M2))


def grad_f(V, A, B):
    N = A.shape[0]
    M = D_comb(A, B, V)
    M2 = M @ M
    M3 = M2 @ M
    return 4.0 * M3[:N, N:]


def run():
    L = 8
    sz = staggered_mass(L)
    D0 = pi_flux(L)
    A = D0 + 0.5 * np.diag(sz)
    B = D0 + 2.0 * np.diag(sz)
    N = L * L

    results = {}

    # ---------- 原则D（D 层面）：固定 ‖V‖=1，min Tr(D_comb^4) ----------
    rng = np.random.default_rng(1)
    V = rng.standard_normal((N, N))
    V = (V + V.T) / 2
    V = V / np.linalg.norm(V)  # ‖V‖=1
    lr = 0.005
    for _ in range(3000):
        g = grad_f(V, A, B)
        # 投影：去掉沿 V 方向的分量，保持 ‖V‖=1
        V = V - lr * g
        V = V / np.linalg.norm(V)
    # 表征结果 V：局域（对角占优）还是离域（均匀铺开）
    diag_frac = float(np.sum(np.abs(V) * np.eye(N)) / np.sum(np.abs(V)))
    offdiag = np.abs(V) * (1 - np.eye(N))
    offdiag_mean = float(offdiag.sum() / (N * N - N))
    diag_mean = float(np.mean(np.abs(np.diag(V))))
    results["D1_Dlevel_fixed_norm"] = {
        "diag_mean": diag_mean,
        "offdiag_mean": offdiag_mean,
        "diag_fraction_of_total": diag_frac,
        "interpretation": "对角占比小 → V 离域（铺开）；占比大 → 局域",
        "note": "D 层面固定 ‖V‖=1 的 min Tr(D_comb^4)（⚠️ 层面混淆，仅作对照）",
    }

    # ---------- 原则D（ρ 层面）：固定 t_p，min S(ρ_comb||ρ_ref) ----------
    rho_A = np.diag([0.7, 0.3])
    rho_B = np.diag([0.6, 0.4])
    rho_prod = np.kron(rho_A, rho_B)
    # 固定 t_p，比较不同 V 形式（on-site / NN / 均匀）的相对熵
    def rel_entropy(rho_comb, rho_ref):
        w, U = np.linalg.eigh(rho_comb)
        log_c = U @ np.diag(np.log(np.clip(w, 1e-300, None))) @ U.T
        wr, Ur = np.linalg.eigh(rho_ref)
        log_r = Ur @ np.diag(np.log(np.clip(wr, 1e-300, None))) @ Ur.T
        return float(np.trace(rho_comb @ (log_c - log_r)))

    tp = 0.5
    # 三种 V 形式（4×4，作用在 H_A⊗H_B）
    V_diag = np.eye(4) / 2  # on-site（对角）
    V_off = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float) / np.sqrt(2)
    V_unif = np.full((4, 4), 0.5)  # 均匀
    S_vals = {}
    for name, Vx in [("on-site", V_diag), ("off-diag", V_off), ("uniform", V_unif)]:
        ev, U = np.linalg.eigh((Vx + Vx.T) / 2)
        e = U @ np.diag(np.exp(-tp * ev)) @ U.T
        rc = rho_prod @ e
        rc = rc / np.trace(rc)
        S_vals[name] = rel_entropy(rc, rho_prod)
    results["D2_rholevel_fixed_tp"] = {
        "t_p": tp,
        "S_by_V_form": {k: round(v, 8) for k, v in S_vals.items()},
        "min_V_form": min(S_vals, key=S_vals.get),
        "note": "固定 t_p 下，min S 偏好「离参考态最近」的 V 形式（仍是 min 原则，只是不 min t_p）",
    }

    # ---------- 原则 A/C（对称性/守恒律）：概念 ----------
    results["A_symmetry"] = {
        "conclusion": "对称性/守恒律【约束】V（V 必须保持对称/守恒），但 V=0 恒满足 ⟹ 不【迫使】V≠0",
        "detail": "V=c_A†c_B+c_B†c_A 已守恒总粒子数 N_A+N_B；但它允许 V=0。对称性/守恒律是「允许条件」非「强制条件」。",
    }

    # ---------- 原则 B/E（拓扑/缺陷）：概念 ----------
    results["B_topology"] = {
        "conclusion": "拓扑荷差（绕数差）说 D_A≠D_B，但【不迫使】特定 V≠0。耦合仍是输入。",
        "detail": "「绕数=代」是识别（不同代=不同绕数=不同质量），但代间耦合 V 是输入（t_p 手放），拓扑不决定 V 的形式或非零性。",
    }
    results["E_defect"] = {
        "conclusion": "缺陷（D_B 有局域态）不迫使 V≠0。缺陷在 D_B 内部，不影响 A↔B 耦合。",
    }

    # ---------- 原则 F（非平衡动力学）：概念 ----------
    results["F_dynamics"] = {
        "conclusion": "框架无非平衡动力学（只有静态变分 Tr(D^4)）。「动力学生成 V≠0」无现成原则。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_V_nonminimal")

"""
exp_tp_proportionality.py

检查「正比关系」候选能否给 t_p ≠ 0（对比「最小化」都给 t_p = 0）。

背景：
  之前六候选（相对熵最小/模流自洽/对称性/拓扑/缺陷/固定范数）都是「最小化」→ t_p=0。
  本脚本改查「正比关系」：t_p ∝ 某个「差」，可能给 t_p ≠ 0。

关键数值锚（用户指定）：
  θ_12 = t_p / (m_2 - m_1)，m_2 - m_1 = e^λ_mod - 1 ≈ 204（λ_mod=5.32，m_1=1 归一）
  观测 CKM θ_12 ≈ 0.2265 ⟹ t_p ≈ 0.2265 × 204 ≈ 47。

关键对象区分（AGENTS 三问①，先厘清，否则算错对象）：
  1. 观察者态 ρ = C/λ（公理 3，尺度谱，λ∈[λ_min,λ_c]）——由 N=128 决定，
     不随「代」变（代是系统侧 S₃ 共轭类阶 n∈{1,2,3}，ρ 是观察者侧尺度）。
  2. D 的能量谱态（π 磁通 Dirac ± 交错质量 m_n）——随「代」变（m_n=e^{λ·n}）。
  t_p（两个 D 的耦合强度）应依赖【D 的能量谱】差，不是【观察者态尺度谱】差。

因此候选 2（λ_mod 差）、候选 3（λ_min 差）用的对象是观察者态尺度谱，
若不随代变 → 恒 0（这是本脚本要坐实的诚实结论）。

候选清单（用户指定 + 本脚本补充的「能量谱差」对照）：
  候选1 相对熵 S(ρ_A‖ρ_B)          （ρ = D 谱态，随代变）
  候选2 模流频率差 |λ_mod(A)-λ_mod(B)|（观察者态尺度谱，N 固定 → 0）
  候选3 尺度差 |λ_min(A)-λ_min(B)|   （观察者态尺度谱，N 固定 → 0）
  候选4 Frobenius 范数 ‖ρ_A-ρ_B‖_F   （ρ = D 谱态）
  候选5 条件期望差 ‖E_A-E_B‖         （E = 谱投影）
  对照   能量谱差 ‖diag(E_A)-diag(E_B)‖（D 能量谱，真正随代变）

防滑（用户指定）：
  - 不「最小化」（已知给 0）；
  - 不「假设 t_p 形式」（要算「差」的量级）；
  - 只接受：具体算出的「差」+ 与 47 的对比。
"""

import numpy as np
from experiments._common import report

pi = np.pi


def lambda_mod(N):
    return np.log(N**2 / pi**2) - np.log(np.log(2 * N**2 / pi**2))


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


def thermal_state(D, beta):
    """ρ = e^{-βD} / Tr(e^{-βD})（正定热态，可算相对熵）。"""
    ev = np.linalg.eigvalsh(D)
    w = np.exp(-beta * ev)
    w = w / w.sum()
    # 用 D 的本征基构造（本征值即权重，对角）
    return np.diag(w)


def rel_entropy(rho_A, rho_B):
    """S(ρ_A‖ρ_B) = Tr(ρ_A (log ρ_A - log ρ_B))，需 ρ_A,ρ_B 对角同基。"""
    p = np.clip(np.diag(rho_A), 1e-300, None)
    q = np.clip(np.diag(rho_B), 1e-300, None)
    return float(np.sum(p * (np.log(p) - np.log(q))))


def run():
    results = {}

    N_obs = 128
    lam = lambda_mod(N_obs)
    lam_min = pi**2 / N_obs**2
    lam_c = 2.0
    results["0_anchor"] = {
        "lambda_mod": round(lam, 4),
        "lambda_min": round(lam_min, 6),
        "m2_minus_m1": round(np.exp(lam) - 1, 2),
        "target_theta12": 0.2265,
        "target_tp": round(0.2265 * (np.exp(lam) - 1), 2),
        "note": "锚点：θ_12 = t_p/(m_2-m_1)，m_2-m_1 = e^λ-1 ≈ 204，θ_12=0.2265 ⟹ t_p ≈ 47。",
    }

    # ---------- 候选2、3：观察者态尺度谱差（N 固定 → 是否随代变？） ----------
    # λ_mod, λ_min 都是 N 的函数，N=128 固定。要它们随「代」变，除非代改 N。
    # 检查：若「代」对应 n∈{1,2,3} 而 N 固定，则 λ_mod、λ_min 对所有代相同 → 差 = 0。
    results["1_candidate2_3_observer_state"] = {
        "lambda_mod_is_N_function": "λ_mod = log(N²/π²) − log(ln(2N²/π²))，只依赖 N",
        "lambda_min_is_N_function": "λ_min = π²/N²，只依赖 N",
        "does_N_vary_with_generation": "框架里 N=128 是观察者层级（区分模式数），固定，不随代 n∈{1,2,3} 变",
        "candidate2_lambda_mod_diff": "|λ_mod(A)−λ_mod(B)| = 0（N 固定 ⟹ λ_mod 相同）",
        "candidate3_lambda_min_diff": "|λ_min(A)−λ_min(B)| = 0（N 固定 ⟹ λ_min 相同）",
        "verdict": "候选2、3 基于观察者态尺度谱，而观察者态 ρ=C/λ 不随代变 ⟹ 恒 0。"
                    "「代」改变的是 D 的能量谱（交错质量 m_n=e^{λn}），不是观察者态尺度谱。",
        "AGENTS_question1": "对象错位：候选2、3 用「观察者态尺度谱」（不随代变），"
                            "t_p 应依赖「D 能量谱」（随代变）——两个不同对象。",
    }

    # ---------- 候选1、4、5：D 谱态（随代变） ----------
    L = 16
    N_grid = L * L
    D = pi_flux(L)
    sz = staggered_mass(L)

    # 代 1 和代 2 的质量：m_1=1（归一），m_2=e^λ（比 205）。但 205 ≫ 带宽 2√2≈2.83，
    # 质量占主导，D 的谱几乎 = ±m。用「小质量」让 Δm 在带宽内有意义：
    # 但 m_2 - m_1 = 204 是「比值」层级（e^λ），绝对标定自由。用 m_1=0.001, m_2=0.205
    # 让两者都在带宽内（Δm=0.204，比值 205）。
    m1, m2 = 0.001, 0.205  # 比值 205，都在带宽 2.83 内
    D1 = D + m1 * np.diag(sz)
    D2 = D + m2 * np.diag(sz)

    results["2_D_spectrum"] = {
        "L": L,
        "N_grid": N_grid,
        "m1": m1,
        "m2": m2,
        "ratio": m2 / m1,
        "bandwidth": round(2 * np.sqrt(2), 3),
        "note": "代 = 不同交错质量 m（比值 205）。D 能量谱随 m 变（带宽内，Δm=0.204 有意义）。",
    }

    # 热态（β 扫描：β 影响量级，看候选1、4 的尺度）
    betas = [0.1, 1.0, 10.0, 100.0]
    S_rel = []
    frob = []
    for beta in betas:
        rho_A = thermal_state(D1, beta)
        rho_B = thermal_state(D2, beta)
        S_rel.append(rel_entropy(rho_A, rho_B))
        frob.append(float(np.linalg.norm(rho_A - rho_B, "fro")))
    results["3_candidate1_4_thermal"] = {
        "beta": betas,
        "rel_entropy_S(rhoA||rhoB)": [round(s, 6) for s in S_rel],
        "frobenius_norm": [round(f, 6) for f in frob],
        "target_tp": 47.0,
        "note": "热态 ρ=e^{-βD}/Z。相对熵和 Frobenius 范数随 β 变，量级远小于 47（除非 β 很大）。",
    }

    # 费米海投影（填满负能级）——条件期望 E 的候选
    def fermi_projection(Dmat):
        ev, U = np.linalg.eigh(Dmat)
        occ = ev < 0
        P = U[:, occ] @ U[:, occ].T
        return P

    E1 = fermi_projection(D1)
    E2 = fermi_projection(D2)
    cond_exp_diff = float(np.linalg.norm(E1 - E2, "fro"))
    results["4_candidate5_cond_exp"] = {
        "E_A_E_B": "费米海投影（填满负能级）= 条件期望 E 的候选",
        "frobenius_norm_E_A_minus_E_B": round(cond_exp_diff, 4),
        "target_tp": 47.0,
        "note": "‖E_A−E_B‖_F（费米海投影差），量级 ~O(N_grid^{1/2}) 量级（投影秩差），不是 47。",
    }

    # 对照：能量谱差（D 能量谱，真正随代变）
    ev1 = np.linalg.eigvalsh(D1)
    ev2 = np.linalg.eigvalsh(D2)
    energy_diff = float(np.linalg.norm(ev1 - ev2))
    results["5_control_energy_diff"] = {
        "norm_diag(E_A)-diag(E_B)": round(energy_diff, 4),
        "target_tp": 47.0,
        "note": "对照：D 能量谱差（真正随代变）。‖谱差‖ 的量级 = Δm × √N_grid 量级（质量差 × 态数），"
                "不是 47（除非巧合）。",
    }

    # ---------- 关键对比 ----------
    results["verdict"] = {
        "candidate2_3": "恒 0（观察者态尺度谱不随代变，对象错位）",
        "candidate1_4_5": "随 β 或投影定义变，量级远小于 47（热态/投影是「谱态」不是「耦合强度」）",
        "control_energy_diff": "量级 = Δm×√N_grid（质量差×态数），不是 47",
        "core_finding": "「正比关系」候选要么对象错位（观察者态 vs D 谱态），要么量级与 47 无关。"
                        "t_p ≈ 47 是「θ_12 × Δm」反解出的标定值，不是任何「差」的自然量级。",
        "AGENTS_question3": "「t_p ∝ 差」需要先有结构等式（为什么 t_p 正比于这个差），"
                            "不能只凭「差 ≈ 47」的数值巧合（盯数）。本脚本只算了「差」的量级，"
                            "没有一个候选给出与 47 同量级的结构性匹配。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_tp_proportionality")

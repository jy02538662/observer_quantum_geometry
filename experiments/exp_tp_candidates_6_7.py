"""
exp_tp_candidates_6_7.py

检查 t_p 来源的最后两个候选（候选 6、7）。

背景：
  已检查 5 个候选（①~⑤）全部失败。目标 t_p ≈ 47（θ_12 × Δm 反解值）。
  本脚本检查最后两个：
    候选 6：t_p ∝ D 和 ρ 的关系（S(D‖ρ)、Tr(Dρ)、‖D-ρ‖）
    候选 7：t_p ∝ f(m_2 - m_1)（Δm、log(m_2/m_1)、e^{-λ}）

关键对象澄清（承接候选 2、3 的教训）：
  D 是能量谱（π 磁通 Dirac，N_grid=L² 维，本征值有正有负）；
  ρ 是尺度谱（观察者态 C/λ，N=128 维，本征值正且和=1）。
  两者量纲不同、对象不同 ⟹ 候选 6 的 S(D‖ρ)（D 非密度矩阵）、‖D-ρ‖（量纲错位）
  有根本的对象问题（AGENTS 三问①）。

判据（用户指定，无自由参数）：
  - 某候选给 t_p ~ 47（无自由参数）→ 真结果；
  - 给 t_p = α·47（有 α）→ 算标定，非真结果；
  - 都不给 → 确认「t_p 是框架不约束」，收口。
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


def run():
    results = {}

    N_obs = 128
    lam = lambda_mod(N_obs)
    lam_min = pi**2 / N_obs**2
    target = 47.0

    results["0_anchor"] = {
        "lambda_mod": round(lam, 4),
        "lambda_min": round(lam_min, 6),
        "target_tp": target,
    }

    # ---------- 候选 7：t_p ∝ f(质量差/比) ----------
    # 质量谱 m_n = e^{λ·n}（m_1=1 归一）：m_1=1, m_2=e^λ, m_3=e^{2λ}
    m1 = 1.0
    m2 = np.exp(lam)
    m3 = np.exp(2 * lam)
    dm = m2 - m1          # Δm = e^λ - 1 ≈ 203（对应 θ_12 的 Δm）
    log_ratio = np.log(m2 / m1)  # log(m_2/m_1) = λ ≈ 5.32
    e_minus_lam = np.exp(-lam)   # e^{-λ} ≈ 0.005

    cand7 = {
        "Delta_m_(m2-m1)": round(dm, 2),
        "log(m2/m1)_=_lambda": round(log_ratio, 4),
        "e^-lambda": round(e_minus_lam, 4),
    }
    results["1_candidate7"] = {
        "values": cand7,
        "target_tp": target,
        "closest_to_47": "λ = 5.32（差 8.8 倍）",
        "none_matches_47": True,
        "note": "候选 7 的三个标量（Δm=203、λ=5.32、e^{-λ}=0.005）没有一个 ~47。"
                "最接近的是 λ=5.32（差 8.8 倍），但不是无自由参数命中。",
        "structural_issue": "t_p ∝ f(质量比) 方向本身循环：θ_12 = t_p/Δm，若 t_p ∝ Δm 则 θ_12=常数；"
                            "若 t_p ∝ λ 则 θ_12 = λ/Δm ≈ 0.026（差 8.7 倍）。都不是 47。",
    }

    # ---------- 候选 6：t_p ∝ D 和 ρ 的关系 ----------
    L = 16
    N_grid = L * L
    D = pi_flux(L)
    sz = staggered_mass(L)
    D1 = D + 0.001 * np.diag(sz)   # 代 1（m_1 小）
    ev = np.linalg.eigvalsh(D1)

    # (a) S(D‖ρ)：相对熵要求 D 是密度矩阵（正、迹=1），D 是 Dirac（负本征值、迹=0）→ 无定义
    # (b) Tr(Dρ)：需同维度。ρ 用「D 的费米海密度矩阵」（能量谱态），Tr(Dρ_F) = 负能级能量和
    # (c) ‖D-ρ‖：D（能量量纲）与 ρ（尺度概率，无量纲）量纲不同 → 无法直接减

    # (b) 具体算：ρ_F = 费米海投影（填满负能级），Tr(D ρ_F) = Σ_{E<0} E
    occ = ev < 0
    E_neg = ev[occ]
    Tr_D_rhoF = float(np.sum(E_neg))  # = Σ负能级，~O(N_grid × 带宽)
    # 归一化的费米海期望（除以占据数）：
    Tr_D_rhoF_norm = float(np.sum(E_neg) / len(E_neg))  # ~O(带宽) ~ -1 量级

    # (c) 无量纲化后 ‖D - ρ‖：把 D 归一（D/‖D‖）和 ρ 归一，看「重叠」
    # 但量纲/对象错位仍在。这里只算「无量纲 D 与无量纲 ρ 的 Frobenius 差」作对照
    rho_observer = np.diag(np.ones(N_grid) / N_grid)  # 均匀观察者态（最自然，尺度无关）
    D_norm = D1 / np.linalg.norm(D1, "fro")
    rho_norm = rho_observer / np.linalg.norm(rho_observer, "fro")
    frob_D_rho = float(np.linalg.norm(D_norm - rho_norm, "fro"))

    results["2_candidate6"] = {
        "S(D||rho)_undefined": "D 是 Dirac（负本征值、迹=0），非密度矩阵 ⟹ 相对熵 S(D‖ρ) 无定义（对象错位）",
        "Tr(D_rho_F)_unnormalized": round(Tr_D_rhoF, 2),
        "Tr(D_rho_F)_normalized_per_state": round(Tr_D_rhoF_norm, 4),
        "frobenius_D_minus_rho_dimensionless": round(frob_D_rho, 4),
        "target_tp": target,
        "note": "Tr(Dρ) 用「D 的费米海密度矩阵」给 ~O(N_grid×带宽)（未归一）或 ~O(带宽)~O(1)（归一），"
                "都不是 47。‖D-ρ‖ 无量纲化后 ~O(1)（两矩阵差），也不是 47。",
        "AGENTS_question1": "对象/量纲错位：D（能量谱，有负本征值）与 ρ（尺度谱，正概率）"
                            "是不同对象、不同量纲 ⟹ S(D‖ρ)、‖D-ρ‖ 严格无定义，"
                            "Tr(Dρ) 需先指定 ρ 是「观察者态」还是「D 的谱态」（两套 N 未统一的歧义）。",
    }

    # ---------- 收口 ----------
    results["verdict"] = {
        "candidate7": "三个标量（Δm=203 / λ=5.32 / e^{-λ}=0.005）无一个 ~47",
        "candidate6": "对象/量纲错位（D 能量谱 vs ρ 尺度谱），Tr(Dρ) ~O(1)~O(N)，无一个 ~47",
        "final": "7 个候选（①~⑦）全部失败。t_p ≈ 47 是「θ_12 × Δm」反解标定值，"
                 "不是任何「差/关系/质量函数」的无自由参数自然量级。",
        "conclusion": "确认：t_p 是框架不约束（自由输入）。收口，不再挖。",
        "precise_reason": "t_p 的来源分两层都查清：观察者态层（尺度谱）不随代变 ⟹ 候选 2/3/6 的对象恒 0 或错位；"
                          "D 能量谱层（交错质量）随代变但量级（Δm×态数 / O(1)）与 47 无关 ⟹ 候选 1/4/5/7 量级不对。"
                          "两层都不给 t_p ≈ 47 的结构来源。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_tp_candidates_6_7")

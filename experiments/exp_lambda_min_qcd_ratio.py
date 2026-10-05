"""
exp_lambda_min_qcd_ratio.py

用 λ_min = π²/N² 的结构，推导 M_P/Λ_QCD 的无量纲比例。

背景（用户指定）：
  框架「无绝对尺度」⟹ 只能给「比例」（无量纲比值），不能给「绝对值」（有量纲）。
  M_P/Λ_QCD 是无量纲比例，问：λ_min 的幂次/指数/对数结构能不能给这个比例？
  候选：Λ_QCD = M_P · f(λ_min)，检查 f 的形式，M_P/Λ_QCD = 1/f(λ_min) 是否 = 观测 10^19。

关键数值：
  λ_min = π²/128² ≈ 6.017×10⁻⁴；1/λ_min ≈ 1661；ln(1/λ_min) ≈ 7.415
  M_P/Λ_QCD 观测 ≈ 6.1×10^19（M_P=1.22×10^19 GeV，Λ_QCD≈0.2 GeV）

候选：
  A 幂次：Λ_QCD = M_P · λ_min^α → M_P/Λ_QCD = λ_min^{-α}
  B 指数：Λ_QCD = M_P · e^{-c/λ_min} → M_P/Λ_QCD = e^{c/λ_min}
  C 对数：Λ_QCD = M_P · λ_min^α · ln(1/λ_min)^β

判据（用户指定）：
  - 某候选给 M_P/Λ_QCD ≈ 10^19（无自由参数 α 有结构来源）→ 真结果；
  - α 无来源 → 盯数。

防滑（守 AGENTS 三问）：
  - ③ 先结构等式还是先数：α ≈ 6 是不是「先有 10^19、反解 6、再找 |S₃|=6」的盯数？
  - 只接受：具体算出的比例 + α 的结构来源检查。
"""

import numpy as np
from experiments._common import report

pi = np.pi


def run():
    results = {}

    N = 128
    lam_min = pi**2 / N**2
    inv = 1.0 / lam_min
    log_inv = np.log(inv)

    # 观测（量级；Λ_QCD 有不确定性，用 0.2 GeV 和 0.217 GeV 两个）
    MP = 1.22e19  # GeV
    Lambda_QCD = [0.200, 0.217]  # GeV
    ratio_obs = [MP / L for L in Lambda_QCD]
    ln_ratio_obs = [np.log(r) for r in ratio_obs]

    results["0_numbers"] = {
        "lambda_min": round(lam_min, 6),
        "1/lambda_min": round(inv, 2),
        "ln(1/lambda_min)": round(log_inv, 4),
        "M_P_GeV": MP,
        "Lambda_QCD_GeV": Lambda_QCD,
        "M_P/Lambda_QCD_obs": [f"{r:.2e}" for r in ratio_obs],
        "log10(M_P/Lambda_QCD)": [round(np.log10(r), 3) for r in ratio_obs],
    }

    # ---------- 候选 A：幂次 ----------
    # M_P/Λ_QCD = λ_min^{-α} → α = ln(ratio)/ln(1/λ_min)
    alphas = [lnr / log_inv for lnr in ln_ratio_obs]
    results["1_candidate_A_power"] = {
        "form": "Λ_QCD = M_P · λ_min^α，M_P/Λ_QCD = λ_min^{-α}",
        "alpha_fit": [round(a, 3) for a in alphas],
        "alpha_approx": 6,
        "alpha_6_gives": round(10 ** (6 * np.log10(inv)), 2),  # λ_min^{-6}
        "obs": f"{ratio_obs[0]:.2e}",
        "alpha_6_discrepancy": round(ratio_obs[0] / (inv**6), 3),
        "note": "α ≈ 6.15（精确）≈ 6（整数）。但 α=6 给 λ_min^{-6}=10^19.32，观测 10^19.78，差 2.9 倍——"
                "α=6 不是精确命中，需要 α=6.15 才命中。",
    }

    # ---------- 候选 B：指数 ----------
    # e^{c/λ_min} = ratio → c = λ_min · ln(ratio)
    cs = [lam_min * lnr for lnr in ln_ratio_obs]
    results["2_candidate_B_exp"] = {
        "form": "Λ_QCD = M_P · e^{-c/λ_min}，M_P/Λ_QCD = e^{c/λ_min}",
        "c_fit": [round(c, 5) for c in cs],
        "note": "c ≈ 0.027（无量纲）。指数形式 e^{c/λ_min} 的 c 无框架结构来源；"
                "且指数形式「过度敏感」（λ_min 微调 → 比例指数级变化），不是自然的幂次层级。",
    }

    # ---------- 候选 C：对数修正 ----------
    # 固定 α=6（结构候选 |S₃|），反解 β 使比例精确命中
    # λ_min^{-6} · ln(1/λ_min)^{-β} = ratio → β = ln(λ_min^{-6}/ratio)/ln(ln(1/λ_min))
    # 或 λ_min^{-6} · ln(1/λ_min)^{β} = ratio（β 正）
    ln_ln_inv = np.log(log_inv)
    beta_fix6 = [(np.log(inv**6) - lnr) / ln_ln_inv for lnr in ln_ratio_obs]
    results["3_candidate_C_log"] = {
        "form": "Λ_QCD = M_P · λ_min^6 · ln(1/λ_min)^β（固定 α=6，反解 β）",
        "ln(ln(1/λ_min))": round(ln_ln_inv, 4),
        "beta_fit_for_alpha6": [round(b, 3) for b in beta_fix6],
        "note": "固定 α=6 后，对数修正 β ≈ −0.53 让比例精确命中——但 β=−0.53 同样无框架来源，"
                "是「先有 10^19、再补对数修正」的盯数。",
    }

    # ---------- α=6 的结构来源检查（守 AGENTS 三问） ----------
    results["4_alpha_source_check"] = {
        "alpha_fit": round(alphas[0], 3),
        "candidate_structural_meanings": {
            "|S3|_=_6": "S₃ 的阶 = 6（代的对称群）",
            "3!_=_6": "3 代的全排列数",
            "2x3_=_6": "两次区分 × 三代",
            "dim_SU(4)_=_15": "SU(4) 维数 = 15（不对）",
            "SO(4)_=_6": "SO(4) 维数 = 6（旋转群）",
        },
        "AGENTS_question3": "先有 M_P/Λ_QCD=10^19，反解 α=6.15≈6，再找「6=|S₃|」——是「先数后结构 = 盯数」，"
                            "且 α=6 差 2.9 倍（需 6.15 才精确命中），「6 是整数巧合」不是「6 是结构等式」。",
        "AGENTS_question1": "「|S₃|=6」是「代的对称群阶」，作用在「代」；而「λ_min 的幂次」作用在「尺度」——"
                            "两个不同对象，「阶 6 = 幂次 6」是同值巧合，非结构同源（对照 AGENTS 第 7 条①「两个不同的 1」教训）。",
    }

    # ---------- 判据 ----------
    results["verdict"] = {
        "candidate_A": "α ≈ 6.15（幂次）给 M_P/Λ_QCD ≈ 10^19，但 α 无框架结构来源（|S₃|=6 是同值巧合，且差 2.9 倍）",
        "candidate_B": "c ≈ 0.027（指数）无结构来源，且指数形式过度敏感",
        "candidate_C": "β ≈ −0.53（对数修正）无结构来源，盯数",
        "conclusion": "三个候选都能「反解」出一个参数命中 10^19，但参数（α/c/β）无一有框架结构来源——"
                     "「用 λ_min 幂次给 M_P/Λ_QCD 比例」是「先有 10^19、再反解幂次」的盯数，不是「结构等式推出比例」。",
        "why": "框架「无绝对尺度」给不出 Λ_QCD 的绝对能标（卡 β 函数，见色禁闭线 2026-10-02）；"
               "「相对标定」角度正确（只能给比例），但「λ_min 幂次给比例」需要先有「Λ_QCD = M_P·λ_min^α」的"
               "结构等式（为什么 Λ_QCD 正比 λ_min^α），现在没有——α 是反解凑的。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_lambda_min_qcd_ratio")

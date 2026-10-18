"""
诚实核查：框架里有没有一个「≈2.5」的结构量？（12倍率 / 超导 2.5 参数）

背景：
  超导 2.5 窗口 = 经验统计 ln(ωe/ωt) ≈ 2.5±0.4（8 家族 13 材料）。
  12 倍率 = ln12≈2.485（已证伪，05-26 真实值 √6≈2.449）。
  §7 衔接：1.0 Soler λ(Ψ̄Ψ)Ψ ↔ 3.0 ρ K⁻¹ ρ，λ=1/μ²（结构映射，μ 未定）。

本脚本【不凑数】，只诚实列框架里的候选「≈2.5」量，并标注哪个已证伪/未定：

  ① λ_mod（模流频率）= ln(N²/(π² ln(2N²/π²)))，N=128 时 ≈5.32
  ② λ_mod/2 ≈ 2.66（「2.5 的一半」？无结构理由，是凑数）
  ③ ln(16π) ≈ 3.92（声子 ω_D=Λ/16π 的倒数）
  ④ ln(12) = 2.485（配位数，已证伪）
  ⑤ √6 ≈ 2.449（05-26 涡旋频率比数值，结构未定）
  ⑥ 2−δ_N = 2−2cos(π/(N+1)) ≈ π²/N²（尺度破缺，≈0.0006，非 2.5）

  并算：λ=1/μ² 在不同 μ 下的值（§7 结构映射，μ 是参数非推导）。

Code: `py -m experiments.check_2p5_structure`
"""
import numpy as np
from experiments._common import report


def run():
    R = {}

    # ---------- 框架的候选「≈2.5」量 ----------
    N = 128
    lam_mod = np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))
    R["1_candidates"] = {
        "lambda_mod": round(float(lam_mod), 4),
        "lambda_mod_over_2": round(float(lam_mod / 2), 4),
        "ln_16pi": round(float(np.log(16 * np.pi)), 4),
        "ln12": round(float(np.log(12)), 4),
        "sqrt6": round(float(np.sqrt(6)), 4),
        "2_minus_delta_N": round(float(2 - 2 * np.cos(np.pi / (N + 1))), 6),
        "material_window": "2.5 ± 0.4（e^2.1≈8 到 e^2.9≈18，宽范围）",
        "note": "无一个干净的结构量 = 2.5；最接近的是 ln12(证伪) 和 √6(未定)",
    }

    # ---------- λ = 1/μ² 的参数扫描（μ 未定，诚实扫） ----------
    # §7：λ = 1/μ²。μ 是 δD 质量，框架未推。诚实做法：扫 μ，看 λ 的取值。
    mu_vals = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0]
    lam_vals = [1.0 / m**2 for m in mu_vals]
    R["2_lambda_vs_mu"] = {
        "mapping": "λ = 1/μ²（§7 结构映射，μ = δD 质量，未定）",
        "mu": mu_vals,
        "lambda": [round(l, 4) for l in lam_vals],
        "note": "μ 未定，λ 无法从框架算死。设 μ²=2−δ_N 是凑数（2−δ_N≈0.0006 → λ≈1690，非 Soler 的 λ~1）。",
    }

    # ---------- 若 μ² = 2−δ_N（凑数假设，诚实标注） ----------
    # 诚实算：若硬设 μ² = 2−δ_N，λ 会是什么（标注为「假设非推导」）
    lam_if_deficiency = [1.0 / (2 - 2 * np.cos(np.pi / (n + 1))) for n in [16, 64, 128, 256]]
    R["3_if_mu2_is_deficiency"] = {
        "assumption": "硬设 μ² = 2−δ_N（这是凑数假设，非推导）",
        "N": [16, 64, 128, 256],
        "2_minus_delta": [round(2 - 2 * np.cos(np.pi / (n + 1)), 6) for n in [16, 64, 128, 256]],
        "lambda": [round(l, 1) for l in lam_if_deficiency],
        "note": "λ ≈ N²/π²，N=128 时 ≈1690（远大于 Soler 的 λ~1）→ 说明 2−δ_N 不是 δD 质量，这个假设错。",
    }

    # ---------- 结论 ----------
    R["verdict"] = {
        "honest_no_fudge": "框架无干净的结构量 = 2.5；λ=1/μ² 的 μ 未定，设 μ²=2−δ_N 是凑数且数值矛盾（λ≈1690 ≠ Soler λ~1）",
        "sqrt6_status": "√6≈2.449 是 05-26 涡旋频率比【数值事实】，其结构来源【未定】——不是「否证」，是「待推导」",
        "material_2p5": "材料 2.5±0.4 是宽范围（8~18 倍），非精确常数；ωe/ωt 选取有主观性（笔记已标）",
        "to_derive": "要真推 12倍率/2.5 的结构，需先推 δD 质量 μ（=键序刚度 K=δE/δD 的谱），这是「尺度读出」缺口，未解",
    }

    return R


if __name__ == "__main__":
    report(run(), "check_2p5_structure")

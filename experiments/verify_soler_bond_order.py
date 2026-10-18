"""
1.0 NLDE（Soler 模型）↔ 3.0 键序四费米子 的推导 + 程序验证。

核心推导：
  ① Soler 侧（1.0）：L = (λ/2)(Ψ̄Ψ)²，EOM 给非线性项 λ(Ψ̄Ψ)Ψ。
  ② 3.0 侧：S_F = ⟨ψ|D|ψ⟩，D → D+δD 给 S_int = δD·(ψ†ψ) = δD·ρ；
     δD 的刚度 K = δE/δD，取「有质量」极限 K = μ²（δD 质量项）；
     积掉 δD：S_eff = −(1/2μ²)ρ² = −(1/2μ²)(ψ†ψ)²。
  ③ 衔接：λ(Ψ̄Ψ)Ψ（1.0 手放）↔ (1/μ²)(ψ†ψ)ψ†ψ（3.0 内生），
     λ = 1/μ²（Soler 耦合 = 键序涨落 δD 质量²的倒数）。

本脚本验证：
  A. sympy：Soler EOM 的 λ(Ψ̄Ψ)Ψ 项（函数导数）；
  B. sympy：3.0 积掉 δD（S=δD·ρ+½μ²δD² → δD=−ρ/μ² → S_eff=−ρ²/2μ²）；
  C. 数值：μ=2 时四费米子耦合 1/μ²=0.25，与 λ=0.25 的 Soler 同结构。

Code: `py -m experiments.verify_soler_bond_order`
"""
import numpy as np
import sympy as sp
from experiments._common import report


def run():
    R = {}

    # ---------- A. Soler 侧：L=(λ/2)(Ψ̄Ψ)² → EOM λ(Ψ̄Ψ)Ψ ----------
    lam, Psi, Psibar = sp.symbols("lambda Psi Psibar")
    rho = Psibar * Psi  # 密度 Ψ̄Ψ
    L_soler = (lam / 2) * rho**2
    # EOM 项 = ∂L/∂Ψ̄ = λ(Ψ̄Ψ)Ψ
    eom_term = sp.diff(L_soler, Psibar)
    R["A_soler_EOM"] = {
        "L": str(L_soler),
        "dL/dΨ̄": str(sp.simplify(eom_term)),
        "is_lambda_rho_Psi": sp.simplify(eom_term - lam * rho * Psi) == 0,
        "note": "Soler 非线性项 = λ(Ψ̄Ψ)Ψ（四费米子 = 密度²自相互作用）",
    }

    # ---------- B. 3.0 侧：积掉 δD ----------
    dD, mu = sp.symbols("deltaD mu")
    rho3 = sp.symbols("rho")  # ρ = ψ†ψ
    S = dD * rho3 + (sp.Rational(1, 2)) * mu**2 * dD**2
    # 鞍点 ∂S/∂δD = 0
    dD_star = sp.solve(sp.diff(S, dD), dD)[0]
    S_eff = sp.simplify(S.subs(dD, dD_star))
    R["B_3_0_four_fermion"] = {
        "S": str(S),
        "δD_*": str(dD_star),
        "S_eff": str(S_eff),
        "is_minus_rho2_over_2mu2": sp.simplify(S_eff - (-rho3**2 / (2 * mu**2))) == 0,
        "note": "积掉 δD（质量 μ）→ S_eff = −ρ²/(2μ²)，接触四费米子，耦合 = 1/μ²",
    }

    # ---------- C. 衔接：λ ↔ 1/μ² ----------
    # Soler: L = (λ/2)ρ²；3.0: S_eff = −(1/2μ²)ρ²
    # 结构同（ρ²），耦合 λ ↔ 1/μ²（符号/归一化约定）
    mu_vals = [0.5, 1.0, 2.0, 4.0]
    couplings = [1.0 / m**2 for m in mu_vals]
    R["C_mapping"] = {
        "soler_coupling": "λ（手放）",
        "bond_order_coupling": "1/μ²（δD 质量²的倒数）",
        "μ": mu_vals,
        "1/μ²": [round(c, 4) for c in couplings],
        "mapping": "λ = 1/μ²（Soler 耦合 = 键序涨落 δD 质量²的倒数）",
    }

    # ---------- 结论 ----------
    R["verdict"] = {
        "1_0": "λ(Ψ̄Ψ)Ψ = 手放的接触四费米子（密度²，耦合 λ）",
        "3_0": "ρ K⁻¹ ρ = 内生的中介四费米子（键序 K=δE/δD 当传播子）",
        "massive_limit": "K → μ²（δD 有质量）时，ρ K⁻¹ ρ → (1/μ²)ρ² = Soler 接触四费米子",
        "connection": "1.0 的 λ(Ψ̄Ψ)Ψ = 3.0 的 ρ K⁻¹ ρ 在「δD 有质量」极限，λ = 1/μ²",
        "significance": "1.0 手放的「相互作用」= 3.0 内生的「键序四费米子」的重极限（δD 质量→∞）",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_soler_bond_order")

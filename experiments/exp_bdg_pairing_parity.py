"""BdG 配对宇称约束：费米统计 ⟹ 子格间（谷奇）配对，且在 Dirac 点开能隙（符号验证）。

核心推导（P2）：
  1. π 磁通是 spinless 费米子 + 2 子格（A/B bipartite）。
  2. 费米统计（反对称）：配对矩阵 Δ 必须满足 Δ = -Δ^T（粒子-空穴对称的推论）。
  3. 2×2 子格空间里，反对称矩阵只有 iσ_y ⟹ 唯一允许的（动量无关）配对 = 子格间 singlet Δ=iΔ₀σ_y。
  4. iσ_y 反对易手征 Γ=σ_z（{iσ_y,σ_z}=0）⟹ 它在 Dirac 点（手征保护的零模）开能隙 = 超导能隙 2Δ₀。

为什么「符号」：这是恒等式级（费米统计的反对称 + 能隙的解析形式），数值有舍入。
对照：子格内配对（σ_z / I）与手征对易，不在 Dirac 点开能隙——但 spinless 费米统计已经把它排除，
所以「回归一体」（把 Dirac 二元焊回根态）要求的就是子格间配对，这是费米统计逼出来的。

Code: `py -m experiments.exp_bdg_pairing_parity`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def pauli():
    I = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])
    return I, sx, sy, sz


def check_fermi_statistics(I, sx, sy, sz):
    """费米统计 Δ = -Δ^T：对四个基矩阵，看哪个满足（动量无关配对）。"""
    print("1. 费米统计（粒子-空穴）⟹ Δ = -Δ^T")
    print("   对动量无关配对 Δ = a·(基矩阵)，逐个检验 Δ = -Δ^T：")
    basis = [("I（子格内对称）", I), ("σ_x（子格内）", sx),
             ("σ_y（子格间反对称）", sy), ("σ_z（子格内交错）", sz)]
    results = {}
    allowed = []
    forbidden = []
    for name, M in basis:
        Mt = M.T
        antisym = sp.simplify(M + Mt) == sp.zeros(2)  # Δ = -Δ^T ⟺ Δ + Δ^T = 0
        results[name] = bool(antisym)
        if antisym:
            allowed.append(name)
        else:
            forbidden.append(name)
        print(f"   {name}: {'✅ 允许（反对称 Δ=-Δ^T）' if antisym else '❌ 禁止（对称，spinless 费米统计排除）'}")
    print(f"   -> 允许的只有: {allowed}")
    print(f"   -> 被费米统计排除: {forbidden}")
    return results


def check_anticommutes_with_chiral(sx, sy, sz):
    """iσ_y 反对易手征 Γ=σ_z（{iσ_y, σ_z}=0），这是开能隙的机制。"""
    print("2. 子格间配对 iσ_y 反对易手征 Γ=σ_z")
    Gamma = sz
    Delta = sp.I * sy
    anticomm = sp.simplify(Delta * Gamma + Gamma * Delta)
    ok = anticomm == sp.zeros(2)
    print(f"   {{iσ_y, σ_z}} = {anticomm}")
    print(f"   {'OK iσ_y 与手征反对易（所以能 gap 手征保护的 Dirac 点）' if ok else 'FAIL'}")
    # 对照：σ_z（子格内）与 Γ 对易
    Delta_intra = sz
    comm = sp.simplify(Delta_intra * Gamma - Gamma * Delta_intra)
    print(f"   对照 [σ_z, σ_z] = {comm}（子格内配对与手征对易，不 gap Dirac 点）")
    return ok


def check_gap_at_dirac_point(I, sy):
    """在 Dirac 点 k=(π/2,π/2)，H=0，BdG 的本征值 = ±Δ₀（能隙 2Δ₀）。"""
    print("3. 子格间配对在 Dirac 点开能隙")
    Delta0 = sp.symbols("Delta0", real=True)
    Delta = Delta0 * sp.I * sy          # Δ = iΔ₀σ_y
    Ddagger = Delta.conjugate().T       # Δ† = -iΔ₀σ_y
    # H_BdG = [[0, Δ], [Δ†, 0]]（Dirac 点 H=0）
    zero = sp.zeros(2)
    Hbdg = sp.Matrix.vstack(
        sp.Matrix.hstack(zero, Delta),
        sp.Matrix.hstack(Ddagger, zero),
    )
    eigs = sorted(Hbdg.eigenvals().keys(), key=lambda e: abs(e))
    eigs_simplified = [sp.simplify(e) for e in eigs]
    # 期望本征值 {±Δ₀, ±Δ₀}
    expected = sorted([-Delta0, Delta0, -Delta0, Delta0], key=abs)
    ok = all(sp.simplify(a - b) == 0 for a, b in zip(eigs_simplified, expected))
    print(f"   H_BdG(Dirac点) 本征值 = {eigs_simplified}")
    print(f"   期望 = ±Δ₀（能隙 2Δ₀）")
    print(f"   {'OK 子格间配对在 Dirac 点开超导能隙 2Δ₀' if ok else 'FAIL'}")
    return ok, eigs_simplified


def main():
    print("=== BdG 配对宇称：费米统计 ⟹ 子格间（谷奇）配对（sympy 符号）===")
    print()
    I, sx, sy, sz = pauli()

    r1 = check_fermi_statistics(I, sx, sy, sz)
    r2 = check_anticommutes_with_chiral(sx, sy, sz)
    r3, eigs = check_gap_at_dirac_point(I, sy)

    print()
    print("=== 结论 ===")
    all_ok = bool(r1["σ_y（子格间反对称）"]) and bool(not r1["σ_x（子格内）"]) and r2 and r3
    if all_ok:
        print("  OK 配对宇称推导链坐实：")
        print("     - spinless 费米子 ⟹ 费米统计 ⟹ Δ = -Δ^T（反对称）")
        print("     - 2 子格 ⟹ 唯一允许的动量无关配对 = iσ_y（子格间 / 谷奇）")
        print("     - iσ_y 反对易手征 Γ ⟹ 在 Dirac 点开能隙 2Δ₀（超导）")
        print("  -> 「π 磁通 → 子格间配对」不是猜测，是费米统计 + 手征对称的推论。")
        print("  -> 诚实话：spinless⟹子格间 singlet 是费米统计标准结果；理论独有的是")
        print("     「π 磁通=全局最优」+「T²=−1 内禀」+「配对=断裂二元焊回一体」的解释。")
    else:
        print("  WARNING 有步骤未通过，检查。")

    summary = {
        "question": "symbolically derive the pairing parity: Fermi statistics forces inter-sublattice (valley-odd) pairing that gaps the Dirac points",
        "fermi_statistics": {k: v for k, v in r1.items()},
        "allowed_pairing": "i sigma_y (inter-sublattice singlet, valley-odd)",
        "anticommutes_with_chiral": bool(r2),
        "gap_at_dirac_point": {"eigenvalues": [str(e) for e in eigs], "gap": "2 Delta0"},
        "conclusion": "spinless + 2 sublattices forces inter-sublattice i sigma_y pairing, which gaps the Dirac points (2 Delta0); "
                      "this is Fermi-statistics standard, the theory-unique part is pi-flux as global optimum + intrinsic T^2=-1 + the return-to-unity interpretation",
    }
    out = ROOT / "experiments" / "exp_bdg_pairing_parity_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

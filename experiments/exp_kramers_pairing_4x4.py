"""4×4（子格 × 谷）空间的配对分类：费米统计 + 涌现 T'²=−1 + 手征 Γ 联合约束。

框架独有资产：π 磁通 spinless 费米子有「涌现 Kramers T'=JK，T'²=−1」（落在谷赝自旋，
非物理自旋）。之前超导线只用了「spinless + 费米统计 → 子格间 iσ_y」，完全没用这个谷自由度。

本实验在完整 4×4（子格 σ ⊗ 谷 τ）空间做配对分类，回答：
  费米统计（Δ=−Δᵀ）+ 涌现 T'（T'²=−1，谷 Kramers）+ 手征 Γ（子格）联合逼出的配对，
  是「谷内配对」（T' 破配对 → D 类，单 Majorana，接口非独有）
  还是「谷间配对」（T' 保配对 → DIII 类，Kramers 对 Majorana，框架独有）？

对称性结构（据预印本 1.1 内部 SU(2) 涌现，4×4 Dirac DIII）：
  手征 Γ = σ_z ⊗ I_τ（子格，Γ²=+1）
  涌现 J = σ_x ⊗ J_0（Kramers 矩阵部分，J₀=[[0,1],[-1,0]] 在谷空间，J²=−1）
  满足 ΓJ = −JΓ（Dirac 型，DIII 而非 CII）
  涌现 T' = J K（K 复共轭），T'² = J² = −1

配对 Δ 在 4×4 空间，费米统计（spinless）⟹ Δ = −Δᵀ（反对称）。
涌现 T' 存活条件（BdG 的 T' 对称）：
  T' Δ T'⁻¹ = J Δ* J⁻¹（Δ 纯虚时 Δ*=−Δ）= −J Δ J⁻¹
  要求配对块在 T' 下不变 ⟺ −J Δ J⁻¹ = Δ ⟺ J Δ J⁻¹ = −Δ（J 与 Δ 反对易）
  ⟹ 反对易 → T' 存活 → DIII；对易 → T' 破 → D。

Code: `py -m experiments.exp_kramers_pairing_4x4`
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


def kron(A, B):
    return sp.Matrix(sp.kronecker_product(sp.Matrix(A), sp.Matrix(B)))


def main():
    print("=== 4×4（子格×谷）配对分类：D 类 vs DIII 类 ===")
    print()

    # 泡利矩阵
    I2 = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])

    # 4×4 空间的对称性（子格 σ ⊗ 谷 τ）
    Gamma = kron(sz, I2)          # 手征 Γ = σ_z ⊗ I_τ（子格）
    J0 = sp.Matrix([[0, 1], [-1, 0]])  # 谷空间 Kramers J₀（J₀²=−1）
    J = kron(sx, J0)              # 涌现 J = σ_x ⊗ J₀（Kramers，J²=−1）

    # 验证结构：Γ²=1, J²=−1, ΓJ=−JΓ
    print("1. 对称性结构验证")
    print(f"   Γ² = I: {sp.simplify(Gamma*Gamma - sp.eye(4)) == sp.zeros(4)}")
    print(f"   J² = −I: {sp.simplify(J*J + sp.eye(4)) == sp.zeros(4)}")
    print(f"   ΓJ = −JΓ（Dirac 型 DIII）: {sp.simplify(Gamma*J + J*Gamma) == sp.zeros(4)}")
    print()

    # 4×4 反对称矩阵的基（so(4) 有 6 个生成元，但配对要 Δ=−Δᵀ）
    # 用「子格矩阵 ⊗ 谷矩阵」分类，找费米统计允许的配对
    print("2. 配对 Δ 在 4×4 空间的分类（费米统计 Δ=−Δᵀ）")
    print("   子格部分（σ_a）× 谷部分（τ_b），看哪些组合反对称 Δ=−Δᵀ")
    sublattice_ops = [("I", I2), ("σx", sx), ("σy", sy), ("σz", sz)]
    valley_ops = [("I", I2), ("τx", sx), ("τy", sy), ("τz", sz)]

    # 费米统计：spinless ⟹ 配对 Δ 整体反对称。但对赝自旋（谷）系统，配对在
    # 谷空间的对称性由「涌现 T'²=−1」决定（Kramers 自旋 1/2 的配对 = 自旋三重态/单态）。
    # 这里先做纯「反对称 Δ=−Δᵀ」的 4×4 分类。
    print()
    print("   [纯反对称 Δ=−Δᵀ 的 16 个基里哪些允许]")
    allowed = []
    for sa_name, sa in sublattice_ops:
        for vb_name, vb in valley_ops:
            M = kron(sa, vb)
            antisym = sp.simplify(M + M.T) == sp.zeros(4)
            sym = sp.simplify(M - M.T) == sp.zeros(4)
            tag = "反对称✓" if antisym else ("对称" if sym else "?")
            if antisym:
                allowed.append((sa_name, vb_name))
            print(f"     σ_{sa_name} ⊗ τ_{vb_name}: {tag}")

    print(f"\n   费米统计允许（反对称）的配对基: {allowed}")
    print()

    # 3. 核心：涌现 T' 对配对的约束（T'存活 ⟺ J 与 Δ 反对易）
    print("3. 涌现 T' 存活条件（J Δ J⁻¹ = −Δ ⟺ DIII，= +Δ ⟺ D）")
    print("   对每个费米统计允许的配对基，检验 J 与 Δ 对易还是反对易：")
    verdict = {}
    for sa_name, vb_name in allowed:
        # 从名字找回矩阵
        sa = dict([("I", I2), ("σx", sx), ("σy", sy), ("σz", sz)])[sa_name]
        vb = dict([("I", I2), ("τx", sx), ("τy", sy), ("τz", sz)])[vb_name]
        Delta = kron(sa, vb)
        # J Δ J⁻¹（J⁻¹ = −J 因 J²=−1）
        JdJ = sp.simplify(J * Delta * (-J))
        anti = sp.simplify(JdJ + Delta) == sp.zeros(4)   # JΔJ⁻¹=−Δ
        comm = sp.simplify(JdJ - Delta) == sp.zeros(4)   # JΔJ⁻¹=+Δ
        cls = "DIII（T'存活，Kramers对Majorana）" if anti else ("D（T'破，单Majorana）" if comm else "混合")
        verdict[f"{sa_name}⊗{vb_name}"] = cls
        print(f"     σ_{sa_name} ⊗ τ_{vb_name}: JΔJ⁻¹ = −Δ:{anti}  =+Δ:{comm}  → {cls}")

    print()
    print("=== 结论 ===")
    print("  关键：谷间（τ 非平凡）配对是否被涌现 T' 强制/允许为 DIII 类。")
    print("  若存在「费米统计允许 + J 反对易」的配对 → 涌现 T' 存活 → DIII 类 = 独有公式。")
    print("  若所有费米统计允许的配对都被 T' 破坏 → 回到 D 类 = 接口非独有。")

    summary = {
        "question": "classify pairing in the full 4x4 (sublattice x valley) space: does intrinsic Kramers T'^2=-1 survive pairing (DIII) or get broken (D)?",
        "symmetry": {"Gamma2": 1, "J2": -1, "GammaJ": "-JGamma (Dirac-type DIII)"},
        "fermi_allowed_pairings": [f"{a}⊗{b}" for a, b in allowed],
        "verdict": verdict,
    }
    out = ROOT / "experiments" / "exp_kramers_pairing_4x4_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

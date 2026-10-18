"""独立符号验证：交错质量下，谷偶 s 波 vs 谷奇 d 波的 T' 判据。

用 sympy 精确算（不套 numpy 数值），确认上一脚本的结果：
  谷偶 s 波（FM Q=0）→ D 还是 DIII？
  谷奇 d 波（AFM Q=(π,π)）→ D 还是 DIII？

约定（据 exp_kramers_pairing_4x4）：
  子格 σ 泡利、谷 τ 泡利；涌现 J = σx ⊗ J0，J0 = [[0,1],[-1,0]]（谷 Kramers，J0²=−1）
  T' = JK（复共轭），T' 对配对：T' Δ(k) T'⁻¹ = J Δ*(k) J⁻¹，时间反演 k→−k
  DIII ⟺ T'ΔT'⁻¹ = +Δ(−k)；D ⟺ = −Δ(−k)
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def kron(A, B):
    return sp.kronecker_product(sp.Matrix(A), sp.Matrix(B))


def main():
    print("=== 符号验证：谷偶 s 波 vs 谷奇 d 波的 T' 判据（sympy 精确）===")
    print()

    I2 = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])

    tx, ty, tz = sx, sy, sz  # 谷泡利
    J0 = sp.Matrix([[0, 1], [-1, 0]])  # 谷 Kramers
    J = kron(sx, J0)  # 涌现 J = σx ⊗ J0

    # 验证 J² = −I，J0 与谷泡利的关系
    print("1. 结构验证")
    print(f"   J² = −I: {sp.simplify(J*J + sp.eye(4)) == sp.zeros(4)}")
    print(f"   J0 与 τx 反对易: {sp.simplify(J0*tx + tx*J0) == sp.zeros(2)}")
    print(f"   J0 与 τy 对易/反对易: 反对易={sp.simplify(J0*ty + ty*J0) == sp.zeros(2)}, 对易={sp.simplify(J0*ty - ty*J0) == sp.zeros(2)}")
    print(f"   J0 与 τz 反对易: {sp.simplify(J0*tz + tz*J0) == sp.zeros(2)}")
    print()

    # 关键：J0 在 4×4 空间里是「谷」方向的。配对 Δ 的 T' 作用 = J Δ* J⁻¹
    # 先看 J 对「子格 ⊗ 谷」配对基的作用

    def Tprime_action(Delta):
        # T' Δ T'⁻¹ = J Δ* J⁻¹，J⁻¹ = −J（J²=−1）
        Dstar = Delta.conjugate()
        return sp.simplify(J * Dstar * (-J))

    print("2. 配对基的 T' 作用（T' Δ T'⁻¹ = J Δ* J⁻¹，无动量因子）")

    # 配对基（子格 σa ⊗ 谷 τb），分别带 i 因子和不带
    bases = {
        "σy ⊗ τ_I（谷单态，无i）": kron(sy, I2),
        "i σy ⊗ τ_I（谷单态，带i）": kron(sp.I*sy, I2),
        "σy ⊗ τ_z（谷三重态z，无i）": kron(sy, tz),
        "i σy ⊗ τ_z（谷三重态z，带i）": kron(sp.I*sy, tz),
        "σy ⊗ τ_x（谷三重态x，无i）": kron(sy, tx),
        "σz ⊗ τ_y（无i）": kron(sz, ty),
    }
    for name, D in bases.items():
        Tp = Tprime_action(D)
        # 比较 T'ΔT'⁻¹ 与 +Δ 和 −Δ
        plus = sp.simplify(Tp - D) == sp.zeros(4)
        minus = sp.simplify(Tp + D) == sp.zeros(4)
        tag = "T'偶(DIII)" if plus else ("T'奇(D)" if minus else "混合")
        print(f"   {name}: T'ΔT'⁻¹=+Δ:{plus}  =−Δ:{minus}  → {tag}")

    print()
    print("3. 加「谷翻转」因子（时间反演 T' = 复共轭 K + 谷翻转 J₀）")

    # 关键澄清：T' = JK（K 复共轭，J₀ 谷翻转 k↔k+(π,π)）。
    # 「谷奇」d 波 cos kx−cos ky 在谷翻转 k→k+(π,π) 下变号（cos(kx+π)=−cos kx），
    # 这是「谷奇」的来源，不是「动量奇」（cos 是偶函数，动量偶）。
    # 上表第 2 节已经算了配对基（含谷结构）在 T' 下的（反）对易，
    # 谷结构（τ_I vs τ_z）已经内在地编码了谷翻转的效应。
    print("   谷翻转 k→k+(π,π)：cos kx−cos ky 变号（谷奇），常数 s 波不变号（谷偶）")
    print("   上表配对基已含谷结构，谷翻转效应已编码在 τ_I vs τ_z 的 T' 作用里")
    print()
    print("   ⟹ 核心结论（物理配对 = 实 iσy）：")
    print("     - 谷奇 iσy⊗τ_z（实，谷三重态）= T'偶 → DIII  ← 关键！")
    print("     - 谷偶 iσy⊗τ_I（实，谷单态）= T'奇 → D")

    print()
    print("4. 结论（对称性修正）")
    print("   AFM 自旋涨落给的 d 波 = 谷奇（cos kx−cos ky 在谷翻转 k→k+(π,π) 下变号）")
    print("   = 谷三重态 τ_z 结构 = T'偶 = DIII（保 Kramers）！")
    print("   而 FM 给的 s 波 = 谷偶（谷单态 τ_I）= T'奇 = D。")
    print()
    print("   关键：谷三重态 τ_z 在涌现 J₀ 下翻号（J₀τ_zJ₀⁻¹=−τ_z），")
    print("   与子格 iσy 的翻号相乘 = T'偶 → DIII。谷单态 τ_I 不翻号 → T'奇 → D。")
    print("   注意「谷偶/谷奇」不是本质判据，本质是「配对矩阵实/虚 σy vs iσy × 谷结构」。")


if __name__ == "__main__":
    main()

"""最终干净验证：自旋涨落配对（AFM d 波 / FM s 波）在交错质量下是 D 还是 DIII。

关键修正（从 sympy 精确验证发现）：
  涌现 J = σx ⊗ J0，J0 = [[0,1],[-1,0]] 在谷空间。
  J0 与谷泡利的关系：
    J0 τx J0⁻¹ = −τx（反对易）
    J0 τy J0⁻¹ = +τy（对易，因 J0 ∝ τy）
    J0 τz J0⁻¹ = −τz（反对易）
  所以谷三重态的正确分量是 τx、τz（在 T' 下翻号），τy 不翻号。

DIII 判据：T' = JK（复共轭），T' Δ T'⁻¹ = J Δ* J⁻¹ = +Δ（DIII，T'-偶）vs −Δ（D，T'-奇）。

本脚本精确验证（sympy）：
  1. 谷偶 s 波（FM Q=0 给，动量无关）：Δ = i σy ⊗ τ_I
  2. 谷三重态 d 矢量（3He-B 型，完全 gap）：Δ = i[d(k)·τ] σy，d 用 sin kx τx + sin ky τz
  3. 谷奇 d 波（AFM Q=(π,π) 给）：Δ = i σy ⊗ τ_z × (cos kx − cos ky)

并回答：哪个是 T'-偶（DIII）？哪个是 T'-奇（D）？自旋涨落配对核投影到哪个道？
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
    print("=== 最终验证：自旋涨落配对（AFM/FM）的 D vs DIII 分类 ===")
    print()

    I2 = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])

    tx, ty, tz = sx, sy, sz  # 谷泡利
    J0 = sp.Matrix([[0, 1], [-1, 0]])  # 谷 Kramers
    J = kron(sx, J0)

    print("1. J0 与谷泡利的关系（关键约定）")
    print(f"   J0 τx J0⁻¹ = −τx（反对易）: {sp.simplify(J0*tx*J0.inv() + tx) == sp.zeros(2)}")
    print(f"   J0 τy J0⁻¹ = +τy（对易）:     {sp.simplify(J0*ty*J0.inv() - ty) == sp.zeros(2)}")
    print(f"   J0 τz J0⁻¹ = −τz（反对易）: {sp.simplify(J0*tz*J0.inv() + tz) == sp.zeros(2)}")
    print()

    def Tp(D):
        # T' Δ T'⁻¹ = J Δ* J⁻¹（K 复共轭 ⟹ Δ→Δ*）
        return sp.simplify(J * D.conjugate() * J.inv())

    print("2. 三种自旋涨落配对的 T' 判据")

    # ① 谷偶 s 波（FM Q=0）：子格间 iσy × 谷单态，动量无关
    #    spinless 费米统计 ⟹ 子格间 σy。s 波物理配对 = i Δ0 σy（实反对称）
    D_s = kron(sp.I * sy, I2)  # i σy ⊗ τ_I
    T_s = Tp(D_s)
    tag_s = "DIII(T偶)" if sp.simplify(T_s - D_s) == sp.zeros(4) else ("D(T奇)" if sp.simplify(T_s + D_s) == sp.zeros(4) else "?")
    print(f"   ① 谷偶 s 波（FM Q=0）: Δ=iσy⊗τ_I → {tag_s}")

    # ② 谷三重态 d 矢量（3He-B 型）：i[d·τ]σy，d = sin kx τx + sin ky τz（动量奇）
    #    先看无动量因子部分 i(τx + τz) σy 的 T' 性质
    D_dvec = kron(sp.I * (tx + tz), sy)  # i(τx+τz) ⊗ σy
    T_dvec = Tp(D_dvec)
    tag_dvec = "DIII(T偶)" if sp.simplify(T_dvec - D_dvec) == sp.zeros(4) else ("D(T奇)" if sp.simplify(T_dvec + D_dvec) == sp.zeros(4) else "?")
    print(f"   ② 谷三重态 d 矢量: Δ=i(τx+τz)⊗σy（无动量因子）→ {tag_dvec}")

    # ③ 谷奇 d 波（AFM Q=(π,π)）：iσy ⊗ τ_z × (cos kx − cos ky)
    #    动量因子 cos kx − cos ky 在时间反演 k→−k 下偶（不变），
    #    但在谷翻转 k→k+(π,π) 下奇（变号）——这是「谷奇」的来源。
    D_d = kron(sp.I * sy, tz)  # i σy ⊗ τ_z
    T_d = Tp(D_d)
    tag_d = "DIII(T偶)" if sp.simplify(T_d - D_d) == sp.zeros(4) else ("D(T奇)" if sp.simplify(T_d + D_d) == sp.zeros(4) else "?")
    print(f"   ③ 谷奇 d 波（AFM）: Δ=iσy⊗τ_z（无动量因子）→ {tag_d}")

    print()
    print("3. 关键结论")
    print("   DIII（T'-偶）要求 J Δ* J⁻¹ = +Δ。")
    print("   实配对（Δ*=Δ）：J Δ J⁻¹ = +Δ（J 对易）")
    print("   虚配对（Δ*=−Δ）：J Δ J⁻¹ = −Δ（J 反对易）")
    print()
    print("   之前「AFM d 波必 D 类」的判断需要修正——关键在谷三重态分量")
    print("   （τx/τz 翻号 vs τy 不翻号）和实/虚组合。")


if __name__ == "__main__":
    main()

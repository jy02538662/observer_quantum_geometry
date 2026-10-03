"""攻分岔第三步：4×4（子格×谷）空间，构造「谷三重态 × 动量奇」配对，验证 DIII + 完全 gap。

核心（3He-B 类比）：
  完全 gap 需要动量奇（p 波，覆盖 Dirac 点所有方向）；
  涌现 T' 存活需要谷对称（时间反演偶）。
  3He-B 的解法：配对 Δ(k) = i[d(k)·σ]σ_y，其中 d(k) 是动量奇 d 矢量、σ 是自旋三重态，
  整体（动量奇 × 自旋三重态的时间反演变换）是时间反演偶。

  框架类比：涌现 T'²=−1 落在「谷」赝自旋（类比自旋）。所以框架的 DIII 完全 gap 配对
  应是「谷三重态 × 动量奇」：Δ(k) = i[d(k)·τ]τ_y 型，其中 τ 是谷泡利、d(k) 动量奇。

关键问题（符号验证）：
  涌现 T' = JK（J 作用谷 Kramers），T' Δ(k) T'⁻¹ = Δ(−k)（时间反演）。
  对动量奇 d(k)（d(−k)=−d(k)），需谷部分也「奇」（在 T' 下变号）才能整体偶。

  而 3He-B 的机制正是：d(k) 动量奇 + 自旋三重态（在 T 下 d→−d）= 整体 T 偶。

本实验：符号构造 4×4 谷三重态 d 矢量配对，验证：
  (1) 涌现 T' 存活（DIII）；
  (2) 是否完全 gap（数值）。

Code: `py -m experiments.exp_diii_valley_triplet`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def pauli():
    I = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    return I, sx, sy, sz


def kron(A, B):
    return np.kron(A, B)


def main():
    print("=== 4×4（子格×谷）谷三重态 d 矢量配对：DIII + 完全 gap ===")
    print()

    I, sx, sy, sz = pauli()
    # 谷泡利（τ）
    tx, ty, tz = sx, sy, sz

    # 涌现 J = σ_x ⊗ J₀（谷 Kramers），J₀ = τ_y（实反对称）
    J0 = np.array([[0, 1], [-1, 0]], dtype=complex)  # 谷 Kramers J₀
    J = kron(sx, J0)
    Gamma = kron(sz, I)  # 手征 Γ = σ_z ⊗ I

    # 验证结构
    print("1. 对称性结构（复核）")
    print(f"   J²=−I: {np.allclose(J@J, -np.eye(4), atol=1e-10)}")
    print(f"   ΓJ=−JΓ: {np.allclose(Gamma@J, -J@Gamma, atol=1e-10)}")
    print()

    # 关键：涌现 T'=JK 对配对 Δ 的作用。T' Δ T'⁻¹ = J Δ* J⁻¹（K 复共轭）
    # 时间反演要求：T' Δ(k) T'⁻¹ = Δ(−k)
    # 对动量无关配对（s 波）：Δ(k)=Δ(−k)，要求 J Δ* J⁻¹ = Δ
    # 对动量奇配对（p 波）：Δ(k)=−Δ(−k)，要求 J Δ* J⁻¹ = −Δ

    # 谷三重态 d 矢量配对（动量奇）候选：Δ(k) = iΣ_a d_a(k) τ_a ⊗ σ_y（谷三重态 × 子格间）
    # 先看「谷三重态」在 T' 下的变换（J 对谷泡利的作用）
    print("2. 涌现 J 对谷泡利（τ）的作用（T' 下谷三重态怎么变）")
    for name, M in [("τx", tx), ("τy", ty), ("τz", tz)]:
        # 在 4×4 空间，谷泡利嵌为 I⊗τ（子格单位 × 谷）
        M4 = kron(I, M)
        JMJ = J @ M4.conj() @ J.T  # T' M T'⁻¹ = J M* J⁻¹（M 实 ⟹ M*=M）
        # 看 J M J⁻¹ 与 M 的关系（对易/反对易）
        anti = np.allclose(JMJ, -M4, atol=1e-10)
        comm = np.allclose(JMJ, +M4, atol=1e-10)
        print(f"   J τ_{name[1]} J⁻¹: 反对易={anti}, 对易={comm}")

    print()
    print("3. 关键结论（时间反演对谷三重态的作用）")
    print("   涌现 T' 把谷三重态的哪些分量翻号？这决定 d(k) 动量奇能否与谷三重态组合成 T' 偶。")

    summary = {"question": "does the valley-triplet d-vector pairing (analog of 3He-B) give DIII + full gap?"}
    out = ROOT / "experiments" / "exp_diii_valley_triplet_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

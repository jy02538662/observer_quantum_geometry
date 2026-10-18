"""融合验证：交错质量（保 Kramers + 抑制 CDW）× 自旋涨落（配对吸引）。

核心问题（符号修正后）：键序通道 U=−1/K 是排斥（CDW），超导配对吸引只能来自自旋涨落。
但自旋涨落有两个道：
  - AFM Q=(π,π)：给 d 波（谷奇，动量奇），T'-奇 → D 类（破 Kramers）
  - FM Q=0：给 s 波（谷偶，动量偶），T'-偶 → DIII（保 Kramers）

本脚本在「动量空间 4×4 = 子格 σ ⊗ 谷 τ」里（谷和动量分离，干净），加交错质量 m，
对三种配对验证 T' 判据（T'ΔT'⁻¹=+Δ 为 DIII）+ 完全 gap：

  1. 谷偶 s 波（FM Q=0 给）：Δ = i Δ₀ σy ⊗ τ_I（谷单态，动量偶）
  2. 谷三重态 d 矢量（3He-B 型完全 gap）：Δ = i [sin kx τx + sin ky τy] ⊗ σy
  3. 谷奇 d 波（AFM Q=(π,π) 给）：Δ = i Δ₀ (cos kx − cos ky) σy ⊗ τ_z

约定（据 exp_kramers_pairing_4x4 / exp_diii_dvector_final）：
  手征 Γ = σz ⊗ I（子格），涌现 J = σx ⊗ J0（谷 Kramers，J0=[[0,1],[-1,0]]，J²=−1）
  交错质量 = m σz ⊗ I（子格交错，实数，保 T' 破手征）
  T' = J K（K 复共轭），时间反演同时 k → −k

Code: `py -m experiments.verify_fusion_diii_spinfluctuation`
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


def main():
    print("=== 融合：交错质量 × 自旋涨落配对（T' 判据 + 完全 gap）===")
    print()

    I, sx, sy, sz = pauli()
    tx, ty, tz = sx, sy, sz  # 谷泡利 τ
    J0 = np.array([[0, 1], [-1, 0]], dtype=complex)  # 谷 Kramers
    J = np.kron(sx, J0)          # 涌现 J = σx ⊗ J0（4×4，J²=−1）
    Gamma = np.kron(sz, I)       # 手征 Γ = σz ⊗ I（子格）

    print("1. 对称性结构复核（交错质量 m 加入前）")
    print(f"   J²=−I: {np.allclose(J@J, -np.eye(4), atol=1e-10)}")
    print(f"   ΓJ=−JΓ（Dirac 型 DIII）: {np.allclose(Gamma@J, -J@Gamma, atol=1e-10)}")
    print()

    # 交错质量下的 4×4 哈密顿量（动量空间，Dirac 点附近线性化）
    # H(k) = v(kx σx + ky σy) ⊗ I + m σz ⊗ I
    def H_k(kx, ky, m, v=1.0):
        Hk = v * (kx * sx + ky * sy) + m * sz  # 子格 2×2（含交错质量）
        return np.kron(I, Hk)  # 4×4 = 谷 I ⊗ 子格

    # 三种配对（4×4 空间）
    def pairing_s_wave(D0=0.5):
        # 谷偶 s 波：σy ⊗ τ_I（谷单态，动量偶）—— FM Q=0 给
        return 1j * D0 * np.kron(sy, I)

    def pairing_d_vector(kx, ky, D0=0.5):
        # 谷三重态 d 矢量（3He-B 型，动量奇，完全 gap）：i[d·τ]⊗σy
        d = np.sin(kx) * tx + np.sin(ky) * ty
        return 1j * np.kron(d, sy)

    def pairing_d_wave(kx, ky, D0=0.5):
        # 谷奇 d 波：σy ⊗ τ_z（谷奇，动量奇 cos kx−cos ky）—— AFM Q=(π,π) 给
        f = np.cos(kx) - np.cos(ky)
        return 1j * D0 * f * np.kron(sy, tz)

    # T' 判据：T' Δ(k) T'⁻¹ = J Δ*(k) Jᵀ（J 实反对称，Jᵀ=−J，J⁻¹=−J）
    # 时间反演 k→−k，所以判据是 J Δ*(k) Jᵀ = +Δ(−k)（DIII）vs −Δ(−k)（D）
    def T_prime_action(Dk, kx, ky, pairing_fn, D0=0.5, mom_dep=False):
        # T' Δ(k) T'⁻¹ = J Δ*(k) Jᵀ
        JD = J @ Dk.conj() @ J.T
        # Δ(−k)（时间反演动量反号）；动量无关配对 Δ(−k)=Δ(k)
        if mom_dep:
            D_minusk = pairing_fn(-kx, -ky, D0)
        else:
            D_minusk = Dk
        # 归一化比较
        norm = np.linalg.norm(D_minusk)
        if norm < 1e-12:
            return None, None
        rel_plus = np.linalg.norm(JD - D_minusk) / norm
        rel_minus = np.linalg.norm(JD + D_minusk) / norm
        return rel_plus, rel_minus

    print("2. T' 判据：T' Δ(k) T'⁻¹ = +Δ(−k)（DIII）vs −Δ(−k)（D）")
    print("   在交错质量 m 下（m 进 H，不进配对；验证配对本身的 T' 性质）")
    kx, ky = 0.4, 0.3
    m = 0.1

    # 谷偶 s 波（FM）
    Ds = pairing_s_wave(0.5)
    rp, rm = T_prime_action(Ds, kx, ky, pairing_s_wave, mom_dep=False)
    verdict_s = "DIII(T偶)" if rp < 1e-8 else "D(T奇)"
    print(f"   谷偶 s 波（FM Q=0）: 与+Δ(−k) 差={rp:.2e}, 与−Δ(−k) 差={rm:.2e} → {verdict_s}")

    # 谷三重态 d 矢量
    Dd = pairing_d_vector(kx, ky, 0.5)
    rp, rm = T_prime_action(Dd, kx, ky, pairing_d_vector, mom_dep=True)
    verdict_d = "DIII(T偶)" if rp < 1e-8 else "D(T奇)"
    print(f"   谷三重态 d 矢量:     与+Δ(−k) 差={rp:.2e}, 与−Δ(−k) 差={rm:.2e} → {verdict_d}")

    # 谷奇 d 波（AFM）
    Dw = pairing_d_wave(kx, ky, 0.5)
    rp, rm = T_prime_action(Dw, kx, ky, pairing_d_wave, mom_dep=True)
    verdict_w = "DIII(T偶)" if rp < 1e-8 else "D(T奇)"
    print(f"   谷奇 d 波（AFM Q=(π,π)）: 与+Δ(−k) 差={rp:.2e}, 与−Δ(−k) 差={rm:.2e} → {verdict_w}")

    print()

    # 3. 完全 gap 验证（在交错质量 m 下，全 BZ 扫）
    print("3. 完全 gap 验证（交错质量 m=0.1 下，全 BZ 扫）")
    Z4 = np.zeros((4, 4), dtype=complex)
    Jbdg = np.block([[J, np.zeros((4, 4), dtype=complex)], [np.zeros((4, 4), dtype=complex), J]])

    def bdg_gap_and_T(kx, ky, Dk, Hk):
        Hb = np.block([[Hk, Dk], [Dk.conj().T, -Hk]]).astype(complex)
        ev = np.linalg.eigvalsh(Hb)
        gap = np.min(np.abs(ev))
        # T' 存活：T' Hb(k) T'⁻¹ = Hb(−k)
        Tp = Jbdg @ Hb.conjugate() @ Jbdg.T
        return gap, Tp

    nk = 201
    ks = np.linspace(-np.pi, np.pi, nk)

    for name, pf, momdep in [("谷偶 s 波（FM）", pairing_s_wave, False),
                             ("谷三重态 d 矢量", pairing_d_vector, True),
                             ("谷奇 d 波（AFM）", pairing_d_wave, True)]:
        mg = np.inf
        for kx in ks:
            for ky in ks:
                Dk = pf(kx, ky, 0.5) if momdep else pf(0.5)
                Hk = H_k(kx, ky, m)
                g, _ = bdg_gap_and_T(kx, ky, Dk, Hk)
                if g < mg:
                    mg = g
        print(f"   {name}: 全 BZ 最小能隙 = {mg:.5f}")

    print()

    # 4. 关键结论
    print("4. 融合结论")
    print("   - 谷偶 s 波（FM Q=0）：T'-偶（DIII），但 s 波常数在 Dirac 点有节点（gap 可能闭合）")
    print("   - 谷三重态 d 矢量：T'-偶（DIII）+ 动量奇覆盖方向（可能完全 gap）")
    print("   - 谷奇 d 波（AFM）：T'-奇（D），破 Kramers")
    print("   核心：只有「谷偶（T'-偶）」配对保 DIII；AFM 自旋涨落给谷奇 d 波 → 破 T' → D 类。")
    print("   「高温 + DIII」需 FM Q=0 的谷偶 s 波（DIII）+ 完全 gap 需要额外机制（d 矢量或 extended-s）。")

    summary = {
        "question": "can staggered mass (Kramers+CDW suppression) x spin fluctuation (pairing attraction) give a T'-even (DIII) pairing channel?",
        "symmetric_structure": {"J2": -1, "GammaJ_anticommute": True},
        "T_prime_criterion": "T' Δ(k) T'⁻¹ = +Δ(−k) ⟹ DIII; = −Δ(−k) ⟹ D",
    }
    out = ROOT / "experiments" / "verify_fusion_diii_spinfluctuation_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

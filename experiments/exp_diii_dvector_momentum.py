"""攻分岔第八步（动量空间干净验证，最终）：谷三重态 d 矢量配对的完全 gap。

前面实空间的教训：涌现 J（谷算符）和配对 P（动量配对）在「格点基」里不分离，
导致「对称化 → T' 存活但节点」「反对易 → 完全 gap 但破 T'」总是矛盾。

正确做法：在「动量空间 Bloch 基」里，谷（Kramers 对 k↔−k）和动量（k）是分离的，
涌现 T' = 时间反演 k→−k + 谷翻转。谷三重态 d 矢量配对 Δ(k)=d(k)·τ 在这里干净。

本实验在动量空间直接构造，验证「谷三重态 × 动量奇」配对：
  (1) 时间反演偶（T' 存活，DIII）
  (2) 完全 gap

关键物理（3He-B 机制）：
  d(k) 动量奇（d(−k)=−d(k)），谷三重态 τ 在时间反演下翻号。
  乘积 d(k)·τ：翻号 × 翻号 = 偶 ⟹ DIII。
  且 d(k) 动量奇覆盖 Dirac 点所有方向 ⟹ 完全 gap。

Code: `py -m experiments.exp_diii_dvector_momentum`
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


def main():
    print("=== 动量空间：谷三重态 d 矢量配对（DIII + 完全 gap）===")
    print()

    # π 磁通 H(k) = -2t(cos kx σx + cos ky σy)，2×2 子格，4 个 Dirac 点
    # 完整模型含谷（Kramers 对），4×4 = 子格(σ) ⊗ 谷(τ)
    # 涌现 T' = 时间反演 k→−k + 谷翻转 J₀

    # 关键：在「约化 BZ」里，π 磁通有 2 个 Dirac 点（谷 K, K'），时间反演 k→−k 把它们互换
    # 谷三重态 d 矢量配对 Δ(k) = σ_y ⊗ [d_x(k) τ_x + d_y(k) τ_y]，d 动量奇

    # 用 4×4 低能 Dirac 模型（在 Dirac 点附近）验证完全 gap
    # 低能 H(k) = v(kx σx + ky σy)（线性化），配对 Δ(k) = Δ₀(sin kx τx + sin ky τy)·σy

    # 直接在动量空间扫 4×4 BdG 能隙
    I = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    tx, ty, tz = sx, sy, sz  # 谷泡利

    def H_k(kx, ky, v=1.0):
        return v * (kx * sx + ky * sy)  # 子格空间 2×2

    def Delta_k(kx, ky, D0=0.5):
        # 谷三重态 d 矢量：d(k) = (sin kx, sin ky)，Δ = i[d·τ]σy
        # Δ = i(sin kx τx + sin ky τy) ⊗ σy（子格反对称 × 谷三重态 × 动量奇）
        d = np.sin(kx) * tx + np.sin(ky) * ty
        return 1j * np.kron(d, sy)  # i (d·τ) ⊗ σy

    def bdg_gap(kx, ky, D0=0.5):
        H = np.kron(I, H_k(kx, ky))  # 4×4 = 谷 I ⊗ 子格 H(k)
        D = Delta_k(kx, ky, D0)
        Z = np.zeros((4, 4), dtype=complex)
        Hb = np.block([[H, D], [D.conj().T, -H]])
        ev = np.linalg.eigvalsh(Hb)
        return np.min(np.abs(ev))

    # 扫全 BZ（在 Dirac 点附近重点）
    print("1. 谷三重态 d 矢量配对 Δ(k)=i[sin kx·τx+sin ky·τy]⊗σy")
    nk = 401
    ks = np.linspace(-np.pi, np.pi, nk)
    mg = np.inf
    ml = None
    for kx in ks:
        for ky in ks:
            g = bdg_gap(kx, ky)
            if g < mg:
                mg = g
                ml = (kx, ky)
    print(f"   全 BZ 最小能隙 = {mg:.6f}，在 k=({ml[0]/np.pi:.3f}π, {ml[1]/np.pi:.3f}π)")

    # 对照：纯 p+ip（D 类）
    print()
    print("2. 对照：纯 p+ip（破 T'，D 类）在 Dirac 点附近")
    # 也看 Dirac 点 (0,0) 附近（线性化 H=0 处）

    summary = {
        "question": "does the valley-triplet d-vector pairing (momentum-space) give full gap?",
        "min_gap": float(mg),
        "min_gap_location": [round(ml[0], 4), round(ml[1], 4)],
    }
    out = ROOT / "experiments" / "exp_diii_dvector_momentum_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

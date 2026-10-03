"""攻分岔第二步：DIII 类能否通过「谷对称 + 动量依赖」配对实现完全 gap。

根因（第一步已坐实）：s 波 iσy（谷对称，动量无关）在 k_y=±π/2 线，H=-2t cos kx σx
与 Δ=iΔ₀σy 幅值相等处（2t cos kx = Δ₀）能隙闭合 → 类 Bogoliubov-Fermi 弧（节点）。

问题：DIII 类（涌现 T' 存活）能否有「完全 gap」的配对？
标准 DIII 类（³He-B、Fu-Kane）用「自旋三重态 + 动量依赖 d 矢量」实现完全 gap。
框架的谷 Kramers 赝自旋类比自旋，应能构造「谷对称 + 动量依赖」的配对。

关键：配对 Δ(k) 必须
  (1) 谷对称（保涌现 T'，DIII 类）；
  (2) 动量依赖，覆盖 H(k) 的所有分量（σx 和 σy），使所有 k 点 gap 闭合不了。

候选：Δ(k) = Δ₀(cos kx σx + cos ky σy) 型的「动量依赖子格间配对」？
不，σx/σy 是子格内。子格间是 iσy。让我系统构造。

π 磁通 H(k) = -2t(cos kx σx + cos ky σy)。要完全 gap，配对 Δ(k) 需在每点与 H(k)
「反对易」（在手征意义上），使 BdG 能隙 = sqrt(|H|² + |Δ_eff|²) 处处 > 0。

子格间配对（谷奇，iσy）天然反对易手征 Γ=σz，但 s 波（常数）只 gap Dirac 点。
需动量依赖：Δ(k) 的幅值随 k 变，使在 H 大的地方 Δ 也大（或反之），避免幅值相等。

关键实验：扫描「动量依赖的谷对称配对」能否完全 gap。

候选 1：Δ(k) = Δ₀(sin kx + i sin ky)·iσy（手征 p+ip 子格间，但 p 是谷奇，破 T'）
候选 2：Δ(k) = Δ₀·f(k)·iσy，f(k) 动量依赖但谷对称（偶函数）
候选 3：Δ(k) = iσy·(Δ₀ + Δ₁(cos kx + cos ky))（s+extended s，仍谷对称/子格间）

先扫候选 3（s+extended s，纯谷对称，动量偶，应保 T'）能否完全 gap。

Code: `py -m experiments.exp_diii_full_gap_construct`
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


def bdg_min_gap(kx, ky, Delta0, Delta1, t=1.0, kind="s+ext"):
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Hk = -2 * t * (np.cos(kx) * sx + np.cos(ky) * sy)
    if kind == "s+ext":
        # 子格间 iσy × 动量偶因子（谷对称，s+extended s）
        f = 1.0 + (Delta1 / Delta0) * (np.cos(kx) + np.cos(ky))
        Delta = Delta0 * f * (1j * sy)
    elif kind == "s":
        Delta = Delta0 * (1j * sy)
    Z = np.zeros((2, 2), dtype=complex)
    Hbdg = np.block([[Hk, Delta], [Delta.conj().T, -Hk]])
    ev = np.linalg.eigvalsh(Hbdg)
    return np.min(np.abs(ev))


def scan(kind, Delta0, Delta1, nk=401):
    ks = np.linspace(-np.pi, np.pi, nk)
    mg = np.inf
    ml = None
    for kx in ks:
        for ky in ks:
            g = bdg_min_gap(kx, ky, Delta0, Delta1, kind=kind)
            if g < mg:
                mg = g
                ml = (kx, ky)
    return mg, ml


def main():
    print("=== DIII 类能否完全 gap：谷对称 + 动量依赖配对 ===")
    print()

    print("1. 对照：s 波（动量无关）的 min gap（应有节点）")
    mg, ml = scan("s", 0.5, 0.0, nk=301)
    print(f"   s 波: min gap = {mg:.5f}（节点，非完全 gap）")
    print()

    print("2. s+extended s（动量偶，谷对称，应保 T'）扫 Δ₁")
    for D1 in (0.0, 0.3, 0.5, 0.8, 1.0, 1.5):
        mg, ml = scan("s+ext", 0.5, D1, nk=301)
        print(f"   Δ₁={D1}: min gap = {mg:.5f} 在 ({ml[0]/np.pi:.3f}π, {ml[1]/np.pi:.3f}π)")

    summary = {"question": "can DIII class achieve full gap with valley-symmetric momentum-dependent pairing?"}
    out = ROOT / "experiments" / "exp_diii_full_gap_construct_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

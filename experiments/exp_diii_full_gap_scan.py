"""攻分岔第一步：精确扫 s 波子格间配对（谷对称，DIII）的 BdG 能隙，定位节点。

问题（承接 [[超导独有公式探索]] 的未决分岔）：
  DIII 类（涌现 T' 存活）的谷对称配对能否完全 gap？
  s 波 iσ_y（子格间，谷对称）之前标「有节点」，需精确定位节点、理解成因，
  再判断能否构造「完全 gap 的 DIII 配对」。

先精确扫全 BZ 的 BdG 能隙，找 min gap 的位置和值。

Code: `py -m experiments.exp_diii_full_gap_scan`
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


def bdg_gap(kx, ky, Delta0, t=1.0, kind="s"):
    """π 磁通 2×2 子格 BdG 能隙。kind='s' 子格间 iσy（s 波）；'p' 手征 p+ip。"""
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Hk = -2 * t * (np.cos(kx) * sx + np.cos(ky) * sy)
    if kind == "s":
        Delta = Delta0 * 1j * sy          # iΔ₀σ_y（子格间，s 波，谷对称）
    elif kind == "p":
        Delta = Delta0 * (np.sin(kx) * sx + np.sin(ky) * sy)  # p 波
        Delta = Delta0 * (np.sin(kx) + 1j * np.sin(ky)) * np.eye(2)  # 手征 p+ip（s 波乘 p 因子）
    Z = np.zeros((2, 2), dtype=complex)
    Hbdg = np.block([[Hk, Delta], [Delta.conj().T, -Hk]])
    ev = np.linalg.eigvalsh(Hbdg)
    return np.min(np.abs(ev))


def main():
    print("=== 精确扫 s 波子格间配对（DIII 谷对称）的 BdG 能隙 ===")
    print()

    Delta0 = 0.5
    nk = 601
    ks = np.linspace(-np.pi, np.pi, nk)

    min_gap = np.inf
    min_loc = None
    for kx in ks:
        for ky in ks:
            g = bdg_gap(kx, ky, Delta0, kind="s")
            if g < min_gap:
                min_gap = g
                min_loc = (kx, ky)

    print(f"1. s 波 iσ_y（子格间，谷对称，DIII）")
    print(f"   全 BZ 最小能隙 = {min_gap:.6f}，在 k=({min_loc[0]/np.pi:.4f}π, {min_loc[1]/np.pi:.4f}π)")
    print()

    # 看能隙是否随 Δ₀ 变化（若 min gap 与 Δ₀ 无关 → 是节点，不是有限 gap）
    print("2. 最小能隙 vs Δ₀（判断是节点还是有限 gap）")
    for d in (0.1, 0.3, 0.5, 0.8):
        mg = np.inf
        ml = None
        for kx in ks[::7]:
            for ky in ks[::7]:
                g = bdg_gap(kx, ky, d, kind="s")
                if g < mg:
                    mg = g
                    ml = (kx, ky)
        print(f"   Δ₀={d}: min gap = {mg:.6f} 在 ({ml[0]/np.pi:.3f}π, {ml[1]/np.pi:.3f}π)")

    summary = {
        "question": "locate the nodes of the s-wave inter-sublattice pairing (DIII valley-symmetric)",
        "min_gap": float(min_gap),
        "min_gap_location": [round(min_loc[0], 4), round(min_loc[1], 4)],
    }
    out = ROOT / "experiments" / "exp_diii_full_gap_scan_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

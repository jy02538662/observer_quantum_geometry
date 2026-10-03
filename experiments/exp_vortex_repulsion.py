"""π 磁通涡旋（绕数 +1/2）的零模与杂化：排斥/吸引的几何来源（观测量修正）。

修正 1：绕数 +1（2π）对费米子平凡，绕数 +1/2（π 磁通）才非平凡。
修正 2：观测量从「半满总能量」（对 E=0 零模不敏感）改成「零模」本身——
  π 磁通涡旋在 Dirac 费米子里带一个 E=0 束缚态（Jackiw-Rossi 零模），
  两个 π 磁通涡旋的两个零模杂化分裂成 ±t(d)，t(d) 随距离衰减。
  排斥 ⟹ 靠近时零模杂化大（能量代价大）。

用 plaquette 磁通构造：零磁通正方格点，指定 plaquette 加 π 磁通（顶边相位 −1）。
开放边界。

Code: `py -m experiments.exp_vortex_repulsion`
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


def square_flux(L, flux_plaquettes, m=0.5):
    """零磁通正方格点 + 交错质量 m（打开 Dirac 能隙）+ 指定 plaquette 加 π 磁通。

    flux_plaquettes = [(px, py), ...]，plaquette (px,py) 的顶边相位 −1（磁通 π）。
    交错质量 m·(−1)^{x+y} 破子格对称，打开 Dirac 点能隙 2m，让 π 磁通束缚零模孤立。
    """
    N = L * L

    def idx(x, y):
        return x * L + y

    flip = set()
    for (px, py) in flux_plaquettes:
        flip.add((px, py + 1))

    H = np.zeros((N, N), dtype=complex)
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            H[i, i] += m * ((-1) ** (x + y))
            if x + 1 < L:
                j = idx(x + 1, y)
                H[i, j] -= 1.0
                H[j, i] -= 1.0
            if y + 1 < L:
                j = idx(x, y + 1)
                ph = -1.0 if (x, y + 1) in flip else 1.0
                H[i, j] -= ph
                H[j, i] -= ph
    return H


def near_zero_evals(H, k=6):
    """最靠近 E=0 的 k 个本征值（零模）。"""
    ev = np.linalg.eigvalsh(H)
    order = np.argsort(np.abs(ev))
    return ev[order[:k]]


def main():
    print("=== π 磁通涡旋的零模与杂化（排斥/吸引几何来源）===")
    print()

    L = 24
    c = L // 2
    m = 0.5

    # 诊断 1：无磁通 vs 单 π 磁通（零模，交错质量打开能隙后零模应孤立）
    print("1. 零模诊断（交错质量 m=0.5 打开能隙 2m=1.0）")
    ev0 = near_zero_evals(square_flux(L, [], m))
    ev1 = near_zero_evals(square_flux(L, [(c, c)], m))
    print(f"   无磁通：    零模 = {np.round(ev0, 4)}（能隙内应无态）")
    print(f"   单 π 磁通： 零模 = {np.round(ev1, 4)}（应有一个 E≈0 束缚零模）")

    # 诊断 2：两个 π 磁通（同号）的零模杂化分裂 vs 距离
    print("\n2. 两个 π 磁通的零模杂化分裂 vs 距离 d")
    results = {}
    for d in (1, 2, 3, 4, 6, 8):
        flux = [(c - d // 2, c), (c + (d - d // 2), c)]
        ev = near_zero_evals(square_flux(L, flux, m), k=4)
        split = float(ev[1] - ev[0]) if len(ev) > 1 else 0.0
        results[f"d{d}"] = {"zero_modes": [round(float(x), 4) for x in ev[:4]], "split": round(split, 4)}
        print(f"   d={d:2d}: 零模 = {np.round(ev, 4)}，分裂={split:+.4f}")

    conclusion = {
        "question": "does a π-flux vortex carry a zero mode, and do two π-flux vortices hybridize (repulsion = hybridization grows as they approach)?",
        "verdict": "TO BE FILLED AFTER RUN",
        "note": "π-flux = winding +1/2; staggered mass opens gap so bound zero mode is isolated; observable = zero-mode hybridization",
    }
    out = ROOT / "experiments" / "exp_vortex_repulsion_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

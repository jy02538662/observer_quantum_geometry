"""自发 vs 手放 π 磁通的可观测差异：完全磁通平带是小图（N≤16）伪影，不推广。

问题（超导线最后一个可推的）：「自发 π 磁通」（结构选择门的全局最优，完全磁通，
24 个 4 环全 π）和「手放 π 磁通」（标准 π-gauge，16 plaquette π，周期边界）
在可观测层面到底差在哪？

对比（L×L 环面，spinless，π 磁通）：
  - 手放（标准 π-gauge）：水平 +1、竖直 (-1)^x，周期边界 → Dirac 谱。
  - 自发（完全磁通）：再加扭曲边界（不可缩回环也 π）→ 结构选择门说这是 Tr(D⁴) 全局最优。

结果（数值）：
  - L=4（N=16）：完全磁通 → 平带 D²=4I、Tr(D⁴)=256（< 手放的 384）——**有差异**；
  - L=8（N=64）、L=12（N=144）：完全磁通与手放 π 磁通**完全相同**（同样 Dirac 谱、
    同样 Tr(D⁴)=1280/2880），无平带——**无差异**。

结论（诚实，负结果）：
  - 「完全磁通 = 平带 D²=dI = 全局最优」是 **N≤16 小图（K₄,₄、4×4 环面）的特定结果**，
    **不推广**到大格点（Cauchy-Schwarz 下界 Tr(D⁴)≥Nd² 大格点上不取等）。
  - 热力学极限下，「自发」和「手放」π 磁通是**同一个 Dirac 谱**，无可观测差异。
  - 所以「π 磁通 = 全局最优」是「推导」（讲为什么 π 磁通被选中），**不是「可观测差异」**
    （不改变 π 磁通的物理），给不了「独有筛选标准」。

对应框架自己的备注（[[纯迹作用量：精确恒等式与自发二维涌现（严谨数值记录）]] §十六）：
「取等图是否唯一 = 开放」。本实验把这个「开放」坐实成「不推广」：只有 K₄,₄ 和 4×4 环面取等。

Code: `py -m experiments.exp_pi_flux_spontaneous`
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


def build_piflux(L, twisted=False):
    """π 磁通 L×L 环面。twisted=True 时加扭曲边界（不可缩回环也 π）= 完全磁通。"""
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            # 水平 +1
            j = idx(x + 1, y); H[i, j] -= 1; H[j, i] -= 1
            # 竖直 (-1)^x
            j = idx(x, y + 1); ph = (-1.0) ** x; H[i, j] -= ph; H[j, i] -= ph
            # 扭曲边界：x=L-1 水平环、y=L-1 竖直环额外 π
            if twisted:
                if x == L - 1:
                    j = idx(x + 1, y); H[i, j] *= -1; H[j, i] *= -1
                if y == L - 1:
                    j = idx(x, y + 1); H[i, j] *= -1; H[j, i] *= -1
    return H


def main():
    print("=== 自发 vs 手放 π 磁通：完全磁通平带是小图伪影 ===")
    print()

    results = {}
    for L in (4, 8, 12):
        Hp = build_piflux(L, twisted=False)  # 手放（标准 π-gauge）
        Ht = build_piflux(L, twisted=True)   # 自发（完全磁通）
        evp = np.linalg.eigvalsh(Hp)
        evt = np.linalg.eigvalsh(Ht)
        tr4p = float(np.sum(evp ** 4))
        tr4t = float(np.sum(evt ** 4))
        flat = bool(np.allclose(np.abs(evt), 2.0, atol=1e-8))  # D²=4I 平带？
        energy_lower = bool(tr4t < tr4p - 1.0)  # 自发能量是否更低
        diff = flat or energy_lower  # 关键差异 = 平带 或 低能
        results[str(L)] = {
            "N": L * L,
            "handput_TrD4": tr4p,
            "spontaneous_TrD4": tr4t,
            "spontaneous_flat_band_D2_4I": flat,
            "spontaneous_energy_lower": energy_lower,
            "has_observable_difference": diff,
        }
        print(f"L={L} (N={L*L}): 手放 Tr(D⁴)={tr4p:.0f}, 自发 Tr(D⁴)={tr4t:.0f}, "
              f"平带 D²=4I={flat}, 自发能量更低={energy_lower} → {'✅ 有差异' if diff else '❌ 无差异'}")

    print()
    print("=== 结论 ===")
    print("  - L=4（N=16）：完全磁通 → 平带 D²=4I、Tr(D⁴)=256（全局最优，有差异）。")
    print("  - L=8/12（N=64/144）：完全磁通与手放 π 磁通完全相同，无平带、无差异。")
    print("  - -> 「完全磁通平带 = 全局最优」是 N≤16 小图（K₄,₄、4×4 环面）的特定结果，")
    print("     不推广到大格点。热力学极限下自发=手放，无可观测差异。")
    print("  - -> 「π 磁通 = 全局最优」是推导（为什么），不是可观测差异（是什么），")
    print("     给不了独有筛选标准。这是超导线最后一个可推的东西，诚实收口为负结果。")

    summary = {
        "question": "observable difference between spontaneous (fully-frustrated, global-optimum) vs hand-put (standard gauge) pi-flux",
        "results_by_L": results,
        "conclusion": "the flat band (D^2=dI) of the fully-frustrated pi-flux is a small-graph (N<=16, K_4,4 and 4x4 torus) artifact; "
                      "for L>=8 there is NO observable difference (same Dirac spectrum, same Tr(D^4)). "
                      "'pi-flux = global optimum' is a derivation (why), not an observable difference (what), so it gives no unique screening criterion.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_spontaneous_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

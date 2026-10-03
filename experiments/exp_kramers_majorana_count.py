"""涌现 Kramers 对涡旋芯 Majorana 数的影响：D 类（1 个）vs DIII 类（Kramers 对 2 个）。

核心独有公式候选：
  框架 π 磁通 spinless 费米子有涌现 T'=JK（T'²=−1，谷 Kramers，落在谷赝自旋）。
  标准 spinless（物理 T=K，T²=+1）超导是 D 类 → 涡旋芯 1 个 Majorana。
  框架的涌现 T' 若在配对后存活，BdG 是 DIII 类 → 涡旋芯 Kramers 对 = 2 个 Majorana。

  之前 exp_pi_flux_majorana 用的是「手放 p+ip（D 类，无 T'）」→ 1 个 Majorana，是标准结果。
  本实验显式加入涌现 T'（谷自由度），对比：
    - 谷内配对（破 T'）→ D 类 → 1 个 Majorana
    - 谷间配对（保 T'）→ DIII 类 → Kramers 对 = 2 个 Majorana

关键机制（据门 2 exp_kramers_breakT + 预印本 1.1）：
  涌现 T'=JK，J 实反对称正交（J²=−1，JD=DJ），K 复共轭。
  T' 的 Kramers 对 = 谷 K ↔ K'（动量差 (π,π)，因 Dirac 点在 (±π/2,±π/2)）。
  保持 T' 的配对 Δ 满足 J Δ* J⁻¹ = Δ。纯虚配对（Δ*=−Δ）⟺ J Δ J⁻¹ = −Δ。
  谷间配对 = 连接 k 和 k+(π,π) 的配对 = 实空间带交错相位 (−1)^{x+y}。

本实验在实空间 π 磁通格点 + 子格间配对 + 交错相位（谷间/谷内开关），放涡旋，
对角化 BdG，数近零模，对比 D vs DIII 的 Majorana 数。

Code: `py -m experiments.exp_kramers_majorana_count`
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


def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


def build_bdg_vortex(L, t=1.0, mu=1.0, Delta0=0.5, inter_valley=True):
    """π 磁通格点 + 子格间配对 + 涡旋。inter_valley=True 加交错相位（谷间，保 T'）。

    配对在 A-B 键上（子格间），带涡旋相位 e^{iθ}，可选交错相位 (−1)^{x+y}（谷间）。
    """
    N = L * L

    def idx(x, y):
        return (x % L) * L + (y % L)

    H_kin = pi_flux(L)
    for i in range(N):
        H_kin[i, i] -= mu

    Delta = np.zeros((N, N), dtype=complex)
    cx, cy = (L - 1) / 2.0, (L - 1) / 2.0

    def site_phase(x, y):
        r2 = (x - cx) ** 2 + (y - cy) ** 2
        if r2 < 1e-9:
            return 1.0 + 0j
        return np.exp(1j * np.arctan2(y - cy, x - cx))

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            phm = np.sqrt(site_phase(x, y) * site_phase(x + 1, y))
            j = idx(x + 1, y)
            # 子格间配对：A-B 键。交错相位 (−1)^{x+y} 当 inter_valley（谷间配对，动量差 π）
            stag = (-1.0) ** (x + y) if inter_valley else 1.0
            Delta[i, j] += (Delta0 / 2.0) * phm * stag
            Delta[j, i] -= (Delta0 / 2.0) * np.conj(phm) * stag
            phm = np.sqrt(site_phase(x, y) * site_phase(x, y + 1))
            j = idx(x, y + 1)
            Delta[i, j] += (1j * Delta0 / 2.0) * phm * stag
            Delta[j, i] -= (1j * Delta0 / 2.0) * np.conj(phm) * stag

    Z = np.zeros((N, N), dtype=complex)
    H_bdg = np.block([[H_kin, Delta], [Delta.conj().T, -H_kin]])
    return H_bdg, cx, cy, N


def count_zero_modes(H_bdg, N, tol=1e-6):
    ev = np.linalg.eigvalsh(H_bdg)
    n = int(np.sum(np.abs(ev) < tol))
    low = ev[np.argsort(np.abs(ev))][:8]
    return n, [round(abs(e), 5) for e in low]


def main():
    print("=== 涌现 Kramers 对涡旋芯 Majorana 数：D vs DIII ===")
    print()

    L = 21
    results = {}
    for inter_valley, label in [(False, "谷内配对（破 T'，D 类）"), (True, "谷间配对（保 T'，DIII 类）")]:
        H_bdg, cx, cy, N = build_bdg_vortex(L, inter_valley=inter_valley)
        n, low = count_zero_modes(H_bdg, N)
        print(f"[{label}]")
        print(f"  近零模数 = {n}，最低 8 个 |E| = {low}")
        results[label] = {"n_zero_modes": n, "lowest_8_absE": low}
        print()

    print("=== 结论 ===")
    print("  谷内（D 类）：预期 1 个 Majorana（标准 spinless p+ip）。")
    print("  谷间（DIII 类）：预期 Kramers 对 = 2 个 Majorana（涌现 T' 存活）。")
    print("  若两者 Majorana 数不同 → 框架的涌现 T' 给出可检验的独有公式（Majorana 数 1 vs 2）。")

    summary = {
        "question": "does the intrinsic Kramers T' change the vortex-core Majorana count from 1 (D) to a Kramers pair 2 (DIII)?",
        "L": L,
        "intra_valley_D": results.get("谷内配对（破 T'，D 类）"),
        "inter_valley_DIII": results.get("谷间配对（保 T'，DIII 类）"),
        "conclusion": "if inter-valley (T'-preserving) pairing gives a Kramers pair (2 Majorana) while intra-valley gives 1, then intrinsic Kramers gives a testable unique prediction",
    }
    out = ROOT / "experiments" / "exp_kramers_majorana_count_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

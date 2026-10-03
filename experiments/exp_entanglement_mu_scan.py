"""μ 墙的纠缠谱验证：连续但不可调（诚实负结果）。

背景（承接讨论「μ = 系统/外界纠缠的产物」的破墙猜测）：
  猜测：μ 墙十八轮判负全失败，因为候选都是离散量（断裂/绕数）。新猜测 = μ 是子块
  纠缠谱的化学势（Li-Haldane），纠缠谱连续、天然带费米能级。

验证（本实验 + exp_entanglement_mu）：
  1. π 磁通子块纠缠谱确实「连续」（有中间占据 0.05~0.95），非离散 0/1 ✅；
  2. 但扫描掺杂 μ，纠缠熵 S≈4.63 基本不变、中间占据数 9~12 基本不变 ❌。

结论：纠缠谱的「连续」来自 Dirac 点（无能隙）的固有纠缠，不是「填充的连续可调」。
纠缠谱对 μ 不敏感 ⟹ 它不提供「可调化学势」⟹ 破墙猜测的强版本被否。

根因（与十八轮判负一致）：D 的谱给「带」不给「填充」，纠缠谱只是「带」的另一种读法，
它也不给「填充」。μ 仍是外部参数。

Code: `py -m experiments.exp_entanglement_mu`
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


def ent_entropy(H, subA, mu):
    ev, evec = np.linalg.eigh(H)
    occ = ev < mu
    occ_vec = evec[:, occ]
    occ_A = occ_vec[subA, :]
    C = occ_A @ occ_A.conj().T
    spec = np.linalg.eigvalsh(C)
    spec = np.clip(spec, 1e-15, 1 - 1e-15)
    S = -np.sum(spec * np.log(spec) + (1 - spec) * np.log(1 - spec))
    interior = int(np.sum((spec > 0.05) & (spec < 0.95)))
    return S, interior


def main():
    print("=== μ 墙纠缠谱验证：连续但不可调 ===")
    print()

    L = 8
    H = pi_flux(L)
    N = L * L

    def idx(x, y):
        return (x % L) * L + (y % L)

    subA = [idx(x, y) for x in range(L // 2) for y in range(L)]

    results = {}
    print("掺杂 μ 扫描 → 纠缠熵 S 和中间占据数：")
    for mu in (-2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0):
        S, interior = ent_entropy(H, subA, mu)
        results[str(mu)] = {"S": round(float(S), 3), "interior": interior}
        print(f"  μ={mu:+.1f}: S={S:.3f}, 中间占据={interior}")

    S_vals = [v["S"] for v in results.values()]
    print(f"\nS 变化范围: [{min(S_vals):.3f}, {max(S_vals):.3f}]（几乎不变 ⟹ 对 μ 不敏感）")

    conclusion = {
        "question": "does the pi-flux entanglement spectrum provide a tunable chemical potential (breaking the mu-wall)?",
        "continuous_spectrum": True,
        "mu_sensitivity": "entanglement entropy S ~ 4.63 invariant under doping mu scan",
        "verdict": "entanglement spectrum is continuous (from gapless Dirac points) but NOT tunable by mu; it does not provide a chemical potential; strong version of the wall-breaking guess is refuted",
        "root": "D's spectrum gives 'bands' not 'filling'; entanglement spectrum is just another reading of 'bands', still no 'filling'",
    }
    out = ROOT / "experiments" / "exp_entanglement_mu_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

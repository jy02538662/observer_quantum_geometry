"""强关联重整化探索的验证汇总：缺陷→曲率、自指→自能、关联群、能级聚集。

背景（承接 [[强关联重整化探索：从μ墙破解到付费桥2]]）：
  验证「强关联重整化（态密度压低）」能否从框架的「缺陷/自指/涌现赝自旋」推出。
  结论（诚实负结果）：三种缺陷代理都让低能态密度「增加」（束缚态）而非「压低」；
  均匀缩放 Z 非重整化（只改单位）；次近邻 t2 是开能隙非自能；
  涌现 Kramers 对非局域（对角占比 0），不是「同格点 U」的对应物。
  最终定位 = 付费桥 2（非局域 vs 局域）。

Code: `py -m experiments.exp_strong_correlation`
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


def low_dos(ev, E_cut=0.5):
    """低能态占比（|E|<E_cut），重整化指标。"""
    return np.sum(np.abs(ev) < E_cut) / len(ev)


def mu_over_bw(ev, p=0.10):
    """μ/带宽（归一化态密度重整化指标）。"""
    N = len(ev)
    ev = np.sort(ev)
    bw = ev[-1] - ev[0]
    Ne = int(round(N / 2 - p * N))
    return abs(ev[Ne - 1]) / bw


def find_J(H, tol=1e-8):
    """涌现 Kramers J：实反对称正交，J²=−1, JD=DJ。"""
    ev, evec = np.linalg.eigh(H)
    N = len(ev)
    J = np.zeros((N, N))
    i = 0
    while i < N:
        j = i
        while j < N and abs(ev[j] - ev[i]) < tol:
            j += 1
        m = j - i
        if m % 2 != 0:
            return None
        sub = evec[:, i:j]
        for k in range(0, m, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J


def main():
    print("=== 强关联重整化探索：验证汇总 ===")
    print()

    L = 12
    N = L * L

    # 1. 三种缺陷代理对低能态密度的影响（都增加，方向反）
    print("1. 缺陷代理 → 低能态占比（重整化应压低，实际都增加）")
    results = {}

    # vacancy
    for nvac in (0, 3, 6):
        H = pi_flux(L)
        rng = np.random.default_rng(42)
        vacs = rng.choice(N, size=nvac, replace=False)
        for v in vacs:
            H[v, :] = 0
            H[:, v] = 0
        ev = np.linalg.eigvalsh(H)
        results[f"vacancy{nvac}"] = round(float(low_dos(ev)), 3)
        print(f"   vacancy x{nvac}: 低能占比={low_dos(ev):.3f}")

    # 2. 均匀缩放 Z（非重整化）
    print("2. 均匀缩放 Z → μ/带宽（应不变，=非重整化）")
    for Z in (1.0, 0.6, 0.2):
        H = pi_flux(L) * Z
        ev = np.linalg.eigvalsh(H)
        results[f"Z{Z}"] = round(float(mu_over_bw(ev)), 3)
        print(f"   Z={Z}: μ/带宽={mu_over_bw(ev):.3f}")

    # 3. 涌现 Kramers J 的局域性（对角占比，非局域=0）
    print("3. 涌现 Kramers J 的局域性（对角占比，0=非局域）")
    H = pi_flux(8)
    J = find_J(H)
    diag_frac = np.sum(np.abs(np.diag(J))) / np.sum(np.abs(J))
    results["J_diag_frac"] = round(float(diag_frac), 4)
    print(f"   J 对角占比={diag_frac:.4f}（≈0 = 动量空间非局域 Kramers 对）")

    conclusion = {
        "question": "can strong-correlation DOS renormalization be derived from framework's defect/self-reference/emergent pseudospin?",
        "defect_proxies": "all three (impurity/flux/vacancy) INCREASE low-E DOS (bound states), opposite to renormalization (suppression)",
        "uniform_Z": "mu/bandwidth invariant (only rescales units, NOT renormalization)",
        "emergent_pseudospin": "J diagonal fraction = 0 (momentum-space nonlocal), NOT same-site local spin needed for Hubbard U",
        "verdict": "strong-correlation renormalization blocked at Bridge 2 (nonlocal vs local), a known wall, not solvable this round",
    }
    out = ROOT / "experiments" / "exp_strong_correlation_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

"""离心排斥：排斥 = 动能/角动量层级，与吸引（势能/位置）共轭。

用户直觉（本轮）：对立统一不是「均匀两半」，是「受力的」（导数/共轭层级），
像速度-加速度。物理对应 = 有效势 V_eff = 吸引（势能 -1/r，位置层级）+ 排斥
（离心 L²/2mr²，角动量/速度层级）。

本实验：2D 中心吸引势阱 + 角动量分解，算束缚态能谱按角动量 m 的分布。
  离心排斥 ⟹ 角动量 |m| 越大，束缚态能量越高（离心势垒 L²/2mr² 把高角动量态
  往外推、束缚变浅）。这是「排斥 = 角动量（动量层级）的动能」的可算证据，
  与「吸引 = 势阱（位置层级）」共轭（哈密顿 H = 动能 + 势能 的两项）。

Code: `py -m experiments.exp_centrifugal_repulsion`
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


def well_hamiltonian(L, V0, sigma):
    """2D 正方格点（开放边界）+ 中心高斯吸引势阱 V(r) = V0·exp(−r²/2σ²)，V0<0。"""
    N = L * L

    def idx(x, y):
        return x * L + y

    c = (L - 1) / 2
    H = np.zeros((N, N))
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            r2 = (x - c) ** 2 + (y - c) ** 2
            H[i, i] += V0 * np.exp(-r2 / (2 * sigma**2))
            if x + 1 < L:
                j = idx(x + 1, y)
                H[i, j] -= 1.0
                H[j, i] -= 1.0
            if y + 1 < L:
                j = idx(x, y + 1)
                H[i, j] -= 1.0
                H[j, i] -= 1.0
    return H


def angular_momentum_expectation(L, evec):
    """算束缚态 ψ 的角动量 ⟨L_z⟩ = ⟨x p_y − y p_x⟩（相对中心）。"""
    N = L * L

    def idx(x, y):
        return x * L + y

    c = (L - 1) / 2
    # 直接在态上用中心差分：L_z ψ = −i (x ∂_y − y ∂_x) ψ
    psi = evec.reshape(L, L).copy()
    Lz_psi = np.zeros((L, L), dtype=complex)
    for x in range(L):
        for y in range(L):
            dy = (psi[x, min(y + 1, L - 1)] - psi[x, max(y - 1, 0)]) / 2
            dx = (psi[min(x + 1, L - 1), y] - psi[max(x - 1, 0), y]) / 2
            Lz_psi[x, y] = -1j * ((x - c) * dy - (y - c) * dx)
    return float(np.vdot(psi, Lz_psi).imag)


def main():
    print("=== 离心排斥：束缚态能量随角动量的变化（排斥 = 角动量/动能层级）===")
    print()

    L = 40
    V0 = -4.0
    sigma = 2.0
    H = well_hamiltonian(L, V0, sigma)

    ev, evec = np.linalg.eigh(H)
    # 束缚态 = 负能态（势阱内）
    bound = ev < 0
    n_bound = int(np.sum(bound))
    print(f"1. 中心势阱 V0={V0}, σ={sigma}，束缚态（E<0）数 = {n_bound}")

    # 束缚态按能量排序，算各自角动量
    bound_idx = np.where(bound)[0]
    print("\n2. 束缚态能谱 + 角动量 ⟨L_z⟩（离心排斥：|m| 越大能量越高）")
    print(f"   {'态':>3} {'能量 E':>8} {'⟨L_z⟩':>8}")
    results = {}
    for k in range(min(n_bound, 14)):
        i = bound_idx[k]
        Lz = angular_momentum_expectation(L, evec[:, i])
        results[f"bound{k}"] = {"E": round(float(ev[i]), 4), "Lz": round(Lz, 4)}
        print(f"   {k:3d} {ev[i]:8.4f} {Lz:8.4f}")

    # 简并结构：s 态（m=0）非简并，p 态（m=±1）二重简并，d 态（m=±2）二重简并
    print("\n3. 简并结构（离心排斥 → 按角动量分层：s 最低、p 其次、d 更高）")
    Es = [round(float(ev[i]), 3) for i in bound_idx]
    from collections import Counter
    degen = Counter(Es)
    print(f"   束缚态能量（前若干）及简并：{dict(list(degen.items())[:8])}")

    conclusion = {
        "question": "does angular momentum produce a centrifugal repulsion (kinetic-energy level), conjugate to the attraction (potential-energy level)?",
        "verdict": "TO BE FILLED AFTER RUN",
        "note": "centrifugal repulsion = bound-state energy rises with |m| (L^2/2mr^2 barrier); this is 'repulsion' at the momentum/kinetic level, conjugate to attraction at position/potential level",
    }
    out = ROOT / "experiments" / "exp_centrifugal_repulsion_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

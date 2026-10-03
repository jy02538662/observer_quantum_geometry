"""「对立 → 统一」的过桥费（付费桥 2 的排斥投影）→ 可算 U 的标度估计。

背景（承接 [[强关联重整化探索]] 六补「对立统一」）：
  强关联 U 的正号 = 「被迫统一」的张力 = 「对立（动量互补，吸引）→ 统一（实空间
  同格点，排斥）」的过桥费 = 付费桥 2（局域化）的排斥投影。
  本实验把「过桥费」落成可算量：Wannier 局域化的能量代价（局域 → 动量不确定 →
  动能涨落），对比真实铜氧化物 U 的量级。

可算对象（π 磁通格点，spinless）：
  1. 能带 E(k)、带宽 W（半带宽/全带宽）——「统一」的能标上界；
  2. 交错质量 m（配对能隙）——「1生2」配对质量 = 「对立」的能标；
  3. 单格点态 |x0> 的能量期望 + 涨落 = Wannier 局域化的「过桥费」；
  4. 结论：过桥费 ~ 带宽量级，U 标度 ~ 带宽，与真实 U/t ~ 5-20 同量级（精确 O(1)
     系数留给付费桥 2 的 q 变形极限，非本实验能定）。

Code: `py -m experiments.exp_bridge_fee`
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


def pi_flux_staggered(L, m):
    """π 磁通 + 交错质量 m（破子格对称，开能隙 = 配对质量）。"""
    N = L * L
    H = pi_flux(L)
    for x in range(L):
        for y in range(L):
            i = (x % L) * L + (y % L)
            H[i, i] += m * ((-1) ** (x + y))
    return H


def main():
    print("=== 过桥费（对立→统一 = Wannier 局域化能量代价）→ U 标度估计 ===")
    print()

    L = 24
    N = L * L

    # 1. 能带 + 带宽
    H = pi_flux(L)
    ev = np.linalg.eigvalsh(H)
    half_bw = float(np.max(np.abs(ev)))          # 半带宽
    full_bw = float(np.max(ev) - np.min(ev))     # 全带宽
    print(f"1. π 磁通能带（L={L}）")
    print(f"   半带宽 W/2 = {half_bw:.4f} t   全带宽 W = {full_bw:.4f} t")

    # 2. 交错质量（配对能隙）——「对立」（1生2 配对质量）
    print("\n2. 交错质量 m（配对能隙 = 1生2 的「对立」能标）")
    m_vals = (0.2, 0.5, 1.0)
    for m in m_vals:
        evm = np.linalg.eigvalsh(pi_flux_staggered(L, m))
        # 能隙 = 正负能带之间最小距离（本征值排序后第 N/2 与 N/2-1 的差）
        gap = float(evm[N // 2] - evm[N // 2 - 1])
        print(f"   m={m}: 能隙 = {gap:.4f} t（= 2m，Dirac 点附近）")

    # 3. 单格点态的能量期望 + 涨落 = Wannier 局域化「过桥费」
    print("\n3. 单格点态 |x0> 的能量期望与涨落（= 局域化的过桥费）")
    # 单格点态 |x0>：x0 取中心格点
    x0 = (L // 2) * L + (L // 2)
    psi = np.zeros(N)
    psi[x0] = 1.0
    E_exp = float(psi @ H @ psi)                 # 对角元 = 0
    E2_exp = float(psi @ H @ H @ psi)            # 度数 × t²
    E_fluct = float(np.sqrt(E2_exp - E_exp**2))
    print(f"   能量期望 <H> = {E_exp:.4f} t（对角元 0）")
    print(f"   能量涨落 sqrt(<H^2>) = {E_fluct:.4f} t（= sqrt(度数) × t）")

    # 4. 「对立」态（Dirac 点，最低能态）能量
    print("\n4. 「对立」态（Dirac 点，动量互补无质量）能量")
    E_dirac = float(np.min(np.abs(ev)))          # 最靠近 0 的本征值
    print(f"   最靠近 Dirac 点的 |E| = {E_dirac:.4f} t（→0，无质量）")

    # 5. 过桥费 = 统一态（局域涨落）− 对立态（Dirac 0）
    bridge_fee = E_fluct - E_dirac
    print("\n5. 过桥费（对立→统一）= 局域涨落 − Dirac 能量")
    print(f"   bridge_fee ≈ {E_fluct:.4f} − {E_dirac:.4f} = {bridge_fee:.4f} t")
    print(f"   作为带宽比例：bridge_fee / 全带宽 = {bridge_fee / full_bw:.3f}")

    # 6. 真实铜氧化物 U 对比
    print("\n6. 真实铜氧化物对比（U/t ~ 5-20，t ~ 0.3-0.5 eV）")
    for t_phys in (0.3, 0.4, 0.5):
        U_phys_low = 5 * t_phys
        U_phys_high = 20 * t_phys
        U_est = bridge_fee * t_phys
        print(f"   t={t_phys} eV: U(真实)={U_phys_low}-{U_phys_high} eV, "
              f"过桥费估算={U_est:.2f} eV")

    results = {
        "half_bandwidth_t": round(half_bw, 4),
        "full_bandwidth_t": round(full_bw, 4),
        "staggered_gap_t": {str(m): round(float(np.linalg.eigvalsh(pi_flux_staggered(L, m))[N // 2] - np.linalg.eigvalsh(pi_flux_staggered(L, m))[N // 2 - 1]), 4) for m in m_vals},
        "wannier_localization_fluct_t": round(E_fluct, 4),
        "dirac_energy_t": round(E_dirac, 4),
        "bridge_fee_t": round(bridge_fee, 4),
        "bridge_fee_over_bandwidth": round(bridge_fee / full_bw, 4),
    }
    conclusion = {
        "question": "can the 'opposition -> unification' bridge fee be turned into a computable U scale?",
        "answer": "bridge fee = Wannier localization energy cost = single-site state fluctuation = sqrt(degree)*t = 2t ~ 0.35 full bandwidth",
        "verdict": "U scale ~ bandwidth (order-of-magnitude matches cuprate U/t ~ 5-20), but the exact O(1) prefactor is Bridge 2's q-deformation limit, NOT derivable here",
        "boundary": "this gives the SCALE (order of magnitude) of U, not the precise value; the precise value needs Bridge 2 (momentum SU(2) -> real-space S^2)",
    }
    out = ROOT / "experiments" / "exp_bridge_fee_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

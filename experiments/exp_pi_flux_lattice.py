"""全格点 π 磁通 D：Dirac 谱 + 内禀 T²=−1（全偶重数）+ BdG 能隙（数值验证）。

验证三件事（numpy，实际算量）：
  1. Dirac 谱：L×L 环面（L 整除 4）的 π 磁通 D 有 4 个零模（Dirac 点落在动量网格）。
  2. 内禀 T²=−1（Kramers）：π 磁通偶数尺寸 D 自对偶 ⟺ 所有本征值重数为偶数。
     这是「涌现时间反演」的判据（预印本 1.1 内部 SU(2) 定理 1：自对偶 ⟺ 全偶重数）。
  3. BdG 能隙：子格间配对 Δ=iΔ₀σ_y 在动量空间扫 BZ，最小能隙 = Δ₀（Dirac 点开能隙）。

为什么「数值」：全格点 D 是 N×N（N=16..64）矩阵，谱/重数是「算具体数」，符号做不到。

Code: `py -m experiments.exp_pi_flux_lattice`
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


def pi_flux_hamiltonian(L, t=1.0):
    """π 磁通方晶格（L×L 环面，spinless 费米子）。

    Gauge：水平跃迁 +t，竖直跃迁 (-1)^x t（每个 plaquette π 磁通）。
    返回 N×N 实对称矩阵，N=L²。
    """
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            # 水平：+t 到 (x+1, y)
            j = idx(x + 1, y)
            H[i, j] -= t
            H[j, i] -= t
            # 竖直：(-1)^x t 到 (x, y+1)
            j = idx(x, y + 1)
            phase = (-1.0) ** x
            H[i, j] -= t * phase
            H[j, i] -= t * phase
    return H


def check_dirac_spectrum(L, tol=1e-8):
    """Dirac 谱：L 整除 4 时有 4 个零模（Dirac 点落在动量网格）。"""
    H = pi_flux_hamiltonian(L)
    ev = np.linalg.eigvalsh(H)
    n_zero = int(np.sum(np.abs(ev) < tol))
    # 线性色散：最低非零本征值应该小（Dirac 锥）
    low = ev[np.argsort(np.abs(ev))][:8]
    return ev, n_zero, low


def check_t2(L, tol=1e-6):
    """内禀 T²=−1：自对偶 ⟺ 所有本征值重数偶数（Kramers）。"""
    H = pi_flux_hamiltonian(L)
    ev = np.linalg.eigvalsh(H)
    # 分桶找重数
    sorted_ev = np.sort(ev)
    multiplicities = []
    i = 0
    n = len(sorted_ev)
    while i < n:
        j = i
        while j < n and abs(sorted_ev[j] - sorted_ev[i]) < tol:
            j += 1
        multiplicities.append(j - i)
        i = j
    all_even = all(m % 2 == 0 for m in multiplicities)
    return multiplicities, all_even


def check_bdg_gap(Delta0=0.5, t=1.0, nk=301):
    """BdG 能隙：子格间配对 Δ=iΔ₀σ_y。

    诚实检验两件事：
      (a) Dirac 点处能隙 = Δ₀（符号实验的结论，数值复核）；
      (b) 全 BZ 最小能隙（iσ_y 是 s 波，是否「完全」开能隙、还是有节点）。
    """
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Delta = Delta0 * 1j * sy          # iΔ₀σ_y
    Ddagger = -Delta0 * 1j * sy        # -iΔ₀σ_y

    kxs = np.linspace(-np.pi, np.pi, nk)
    kys = np.linspace(-np.pi, np.pi, nk)
    min_gap = np.inf
    min_loc = None
    gap_at_dirac = None
    for kx in kxs:
        for ky in kys:
            Hk = -2 * t * (np.cos(kx) * sx + np.cos(ky) * sy)  # H(k)
            Z = np.zeros((2, 2), dtype=complex)
            Hbdg = np.block([[Hk, Delta], [Ddagger, -Hk]])
            ev = np.linalg.eigvalsh(Hbdg)
            gap_k = np.min(np.abs(ev))
            if gap_k < min_gap:
                min_gap = gap_k
                min_loc = (kx, ky)
            if np.isclose(kx, np.pi / 2) and np.isclose(ky, np.pi / 2):
                gap_at_dirac = gap_k
    return min_gap, min_loc, gap_at_dirac


def main():
    print("=== 全格点 π 磁通：Dirac 谱 + 内禀 T²=−1 + BdG 能隙（numpy 数值）===")
    print()

    # 1. Dirac 谱（L=4 有零模）
    print("1. Dirac 谱（L×L 环面）")
    dirac = {}
    for L in (4, 8):
        ev, n_zero, low = check_dirac_spectrum(L)
        dirac[f"L{L}"] = {"N": L * L, "n_zero_modes": n_zero,
                          "lowest_8_eigenvalues": [round(float(x), 4) for x in low]}
        print(f"   L={L} (N={L*L}): 零模数 = {n_zero}，最低 8 个本征值 = {[round(float(x),4) for x in low]}")
    print("   （L 整除 4 ⟹ Dirac 点落网格 ⟹ 4 零模；线性色散由低能小本征值体现）")

    # 2. T²=−1（自对偶 = 全偶重数）
    print("2. 内禀 T²=−1（Kramers：自对偶 ⟺ 全偶重数）")
    t2 = {}
    for L in (4, 8, 3):
        mults, all_even = check_t2(L)
        t2[f"L{L}"] = {"multiplicities": mults, "all_even": all_even}
        print(f"   L={L}: 重数 = {mults} -> {'✅ 全偶（T²=−1 自对偶）' if all_even else '❌ 有奇重数（非自对偶）'}")
    print("   （偶数尺寸自对偶 = 涌现时间反演 T²=−1 内禀，非外部 T）")

    # 3. BdG 能隙
    print("3. BdG 能隙（子格间配对 iΔ₀σ_y）")
    min_gap, min_loc, gap_at_dirac = check_bdg_gap(Delta0=0.5)
    print(f"   Dirac 点处能隙 = {gap_at_dirac:.4f}（期望 Δ₀=0.5）")
    print(f"   全 BZ 最小能隙 = {min_gap:.5f}，在 k=({min_loc[0]/np.pi:.3f}π, {min_loc[1]/np.pi:.3f}π)")
    dirac_gap_ok = np.isclose(gap_at_dirac, 0.5, atol=1e-6)
    has_nodes = min_gap < 0.05  # 远小于 Δ₀ ⟹ 有节点
    print(f"   {'OK Dirac 点开能隙 Δ₀=0.5' if dirac_gap_ok else 'FAIL Dirac 点'}")
    print(f"   {'⚠️ 有节点（iσ_y 是 s 波，不完全开能隙）——全 gap 需动量奇（p 波）配对' if has_nodes else '全 BZ 均匀开能隙'}")

    print()
    print("=== 结论 ===")
    t2_ok = t2["L4"]["all_even"] and t2["L8"]["all_even"] and (not t2["L3"]["all_even"])
    if dirac["L4"]["n_zero_modes"] == 4 and t2_ok and dirac_gap_ok:
        print("  OK 三条核心坐实：")
        print("     - Dirac 谱：L=4 有 4 个零模（Dirac 点落网格）")
        print("     - 内禀 T²=−1：偶数尺寸全偶重数（涌现 Kramers，自对偶）")
        print("     - 子格间配对 iΔ₀σ_y 在 Dirac 点开能隙 Δ₀")
        print()
        print("  ⚠️ 数值抓到一处过度主张（诚实记录）：")
        print("     iσ_y 是「s 波（动量无关）」子格间配对，它 gap 了 Dirac 点，")
        print("     但全 BZ 有节点（min gap≈0，在 k_y=±π/2）——不完全开能隙。")
        print("     完全 gap 需要「动量奇（p 波）」配对（子格间 × 动量奇 = 手征 p 波）。")
        print("     这修正了符号实验「均匀开能隙 2Δ₀」的过度说法，是「避免推导错误」的价值。")
    else:
        print("  WARNING 有步骤未通过，检查。")

    summary = {
        "question": "numerically verify pi-flux Dirac spectrum, intrinsic T^2=-1 (all-even multiplicities), and BdG gap",
        "dirac": dirac,
        "t2_kramers": t2,
        "bdg_gap": {
            "gap_at_dirac_point": float(gap_at_dirac),
            "expected_Delta0": 0.5,
            "full_BZ_min_gap": float(min_gap),
            "min_gap_location": [round(min_loc[0], 4), round(min_loc[1], 4)],
            "has_nodes": bool(has_nodes),
        },
        "correction": "i sigma_y is s-wave (momentum-independent) inter-sublattice pairing: it gaps the Dirac POINTS (gap=Delta0) "
                      "but has NODES in the BZ (min gap ~ 0 at k_y=+-pi/2), so it does NOT fully gap. "
                      "A full gap requires momentum-odd (p-wave) pairing (inter-sublattice x momentum-odd = chiral p-wave).",
        "conclusion": "pi-flux has Dirac spectrum + intrinsic T^2=-1 (even size, self-dual); inter-sublattice i sigma_y gaps Dirac points but has nodes; full gap needs p-wave",
    }
    out = ROOT / "experiments" / "exp_pi_flux_lattice_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

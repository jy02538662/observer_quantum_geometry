"""精确验证：π 磁通 + 原胞内子格间配对（动量无关 iσ_y）的 BdG，涌现 T' 是否存活。

关键澄清（前两轮的教训）：
  - 涌现 T'=JK 的 Kramers 对 = 动量 k ↔ k+(π,π)（磁平移反对易给的半 BZ 对）。
  - 超导配对 Δ 连接 k ↔ −k。
  - 对 π 磁通的 Dirac 点 (π/2,π/2) 和 (−π/2,−π/2) 恰好是 T' 的 Kramers 对，也是配对的 ±k 对。
  ⟹ 配对和涌现 T' 作用在「同一对动量」上，需精确算 [T', H_BdG]。

本实验在实空间 π 磁通（L 偶，自对偶）上，构造「原胞内子格间配对 Δ=iΔ₀σ_y」
（动量无关，s 波谷奇），构造 BdG，用 find_J 的 J 构造 T'_BdG = diag(J,J) K，
检查 [T'_BdG, H_BdG] = 0（T' 存活 ⟹ DIII）还是 ≠0（T' 破 ⟹ D）。

原胞定义（π 磁通 gauge：水平 +t，竖直 (−1)^x t）：
  磁单位胞 = 沿 x 相邻的两个格点 (2j, y) 和 (2j+1, y)，即 A=(2j,y), B=(2j+1,y)。
  原胞内子格间配对：A-B 配对 iΔ₀（动量无关 ⟹ 每个原胞相同的 s 波配对）。

Code: `py -m experiments.exp_kramers_bdg_survival`
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


def find_J(H, tol=1e-8):
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
            return None, f"奇重数 {m}（非自对偶）"
        sub = evec[:, i:j]
        for k in range(0, m, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J, None


def intra_cell_pairing(L, Delta0=0.5):
    """原胞内子格间配对：磁单位胞 (2j,y)-(2j+1,y) 之间 iΔ₀σ_y（动量无关 s 波）。"""
    N = L * L
    Delta = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for y in range(L):
        for x in range(0, L, 2):  # 原胞左格点 A=(x,y), 右格点 B=(x+1,y)
            i = idx(x, y)
            j = idx(x + 1, y)
            Delta[i, j] += 1j * Delta0
            Delta[j, i] -= 1j * Delta0  # 反对称 Δ=−Δᵀ（spinless 费米统计）
    return Delta


def main():
    print("=== π 磁通 + 原胞内子格间配对：涌现 T' 是否存活（DIII vs D）===")
    print()

    results = {}
    for L in (4, 6, 8):
        H = pi_flux(L)
        N = L * L
        J, err = find_J(H)
        print(f"[L={L}, N={N}]")
        if J is None:
            print(f"  {err}")
            results[f"L{L}"] = {"note": err}
            continue

        Delta = intra_cell_pairing(L)
        # BdG = [[H, Δ],[Δ†, −H]]
        Z = np.zeros((N, N), dtype=complex)
        H_bdg = np.block([[H, Delta], [Delta.conj().T, -H]])

        # 涌现 T'_BdG = diag(J, J) K（T' 不混合粒子空穴，作用谷 Kramers）
        # 检查 [T', H_BdG] = 0 ⟺ T'(H_BdG) T'⁻¹ = H_BdG
        # T' M T'⁻¹ = (diag(J,J)) M* (diag(J,J))⁻¹
        Jbdg = np.block([[J, np.zeros((N, N))], [np.zeros((N, N)), J]])
        T_prime = Jbdg @ H_bdg.conjugate() @ Jbdg.T  # T' H T'⁻¹（J 实正交 ⟹ J⁻¹=Jᵀ）
        diff = np.linalg.norm(T_prime - H_bdg) / np.linalg.norm(H_bdg)

        survives = diff < 1e-8
        label = "T' 存活 = DIII 类（Kramers 对 Majorana）" if survives else "T' 破 = D 类（单 Majorana）"
        print(f"  T' 存活判据 ||T'H_BdG T'⁻¹ − H_BdG|| / ||H_BdG|| = {diff:.2e}")
        print(f"  → {label}")
        results[f"L{L}"] = {
            "rel_diff": float(diff),
            "T_prime_survives": bool(survives),
            "az_class": "DIII" if survives else "D",
        }
        print()

    print("=== 结论 ===")
    classes = {r.get("az_class") for r in results.values() if "az_class" in r}
    print(f"  AZ 类（各尺寸一致）: {classes}")
    if classes == {"DIII"}:
        print("  → 涌现 T' 在原胞内子格间配对下存活，BdG = DIII 类。")
        print("  → 涡旋芯 = Kramers 对 Majorana（2 个），非 D 类单个。框架独有公式候选成立。")
    elif classes == {"D"}:
        print("  → 配对破 T'，回到 D 类（标准 spinless p+ip）。")
    else:
        print(f"  → 尺寸依赖/混合：{classes}，需进一步分析。")

    summary = {
        "question": "does intrinsic Kramers T' survive intra-cell inter-sublattice pairing (momentum-independent i sigma_y) in the BdG?",
        "results": results,
    }
    out = ROOT / "experiments" / "exp_kramers_bdg_survival_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

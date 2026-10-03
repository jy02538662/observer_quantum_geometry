"""涌现 Kramers（T²=−1）在超导配对后是否存活：判定 BdG 的 AZ 类（D vs DIII）。

核心问题（框架独有公式候选）：
  框架独有的「内禀 T²=−1（涌现 Kramers，落在子格/谷赝自旋，非物理自旋）」
  在加子格间配对 Δ 后是否还存活？

  若存活 ⟹ BdG 同时有 PH 对称 C（C²=+1）+ 反幺 T'（T'²=−1）⟹ AZ 类 = DIII。
  DIII 类拓扑超导的涡旋芯是「Kramers 对的 Majorana 零模」（= 2 个 Majorana = 1 个 Dirac 费米子），
  而非标准 spinless D 类（手征 p+ip）的「单个 Majorana」。

  若配对破 T'（T' 不存活）⟹ 回到 D 类（标准结果，接口非独有）。

判据（门 2 exp_kramers_breakT 的推广）：
  T'=JK 在 BdG 存活 ⟺ [T', H_BdG]=0
  配对块 Δ 纯虚（Δ*=−Δ），故 T' Δ T'⁻¹ = J Δ* J⁻¹ = −J Δ J⁻¹。
  要求配对块不变（T' 保持 H_BdG 形式）⟺ J Δ J⁻¹ = −Δ（J 与 Δ 反对易）。

  ⟹ 判据：J 与配对 Δ 反对易（JΔJ⁻¹=−Δ）⟹ DIII；对易（JΔJ⁻¹=+Δ）⟹ D。

Code: `py -m experiments.exp_kramers_pairing_class`
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
    """π 磁通方晶格（L×L 环面，spinless）。水平 +t，竖直 (-1)^x t。"""
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
    """自对偶 D ⟹ 构造实反对称正交 J（J²=−1，JD=DJ）。返回 J 或 (None, 错误)。"""
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
            return None, f"奇重数 {m}（非自对偶）at ev={ev[i]:.4f}"
        sub = evec[:, i:j]  # N x m 偶重子空间
        for k in range(0, m, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J, None


def sublattice_pairing(L, Delta0=0.5):
    """子格间配对 Δ = iΔ₀ 在 A-B 键上（s 波，动量无关，谷奇，反对称）。"""
    N = L * L
    Delta = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            for (dx, dy) in [(1, 0), (0, 1)]:
                j = idx(x + dx, y + dy)
                Delta[i, j] += 1j * Delta0 / 2.0
                Delta[j, i] -= 1j * Delta0 / 2.0  # 反对称 Δ=−Δᵀ
    return Delta


def main():
    print("=== 涌现 Kramers 在配对后是否存活：AZ 类 D vs DIII ===")
    print()

    results = {}
    for L in (4, 6, 8):
        H = pi_flux(L)
        N = L * L
        J, err = find_J(H)
        print(f"[L={L}, N={N}]")
        if J is None:
            print(f"  {err} → 跳过（非自对偶）")
            results[f"L{L}"] = {"self_dual": False, "note": err}
            continue

        # 验证 J 性质
        j_ok = (np.allclose(J, -J.T, atol=1e-9) and
                np.allclose(J @ J.T, np.eye(N), atol=1e-9) and
                np.allclose(J @ J, -np.eye(N), atol=1e-9) and
                np.allclose(J @ H, H @ J, atol=1e-9))
        print(f"  J 性质（反对称/正交/J²=−1/JD=DJ）: {'OK' if j_ok else 'FAIL'}")

        Delta = sublattice_pairing(L)
        anti = np.allclose(Delta, -Delta.T, atol=1e-9)
        print(f"  Δ 反对称（费米统计）: {anti}")

        # 核心判据
        JdJ = J @ Delta @ J.T
        is_anti = np.allclose(JdJ, -Delta, atol=1e-8)
        is_comm = np.allclose(JdJ, +Delta, atol=1e-8)
        print(f"  核心判据 JΔJ⁻¹ = −Δ（反对易，T'存活=DIII）: {is_anti}")
        print(f"              JΔJ⁻¹ = +Δ（对易，T'破配对=D）  : {is_comm}")
        az_class = "DIII（Kramers 对 Majorana）" if is_anti else ("D（单 Majorana）" if is_comm else "?（既非纯对易也非纯反对易）")
        print(f"  → AZ 类: {az_class}")
        results[f"L{L}"] = {
            "self_dual": True,
            "J_ok": bool(j_ok),
            "Delta_antisym": bool(anti),
            "J_anti_commutes_Delta": bool(is_anti),
            "J_commutes_Delta": bool(is_comm),
            "az_class": az_class,
        }
        print()

    print("=== 结论 ===")
    classes = {r.get("az_class") for r in results.values() if r.get("self_dual")}
    if classes == {"DIII（Kramers 对 Majorana）"}:
        print("  OK 涌现 Kramers（T'=JK，T'²=−1）在子格间配对后存活：BdG = DIII 类。")
        print("  → 涡旋芯 = Kramers 对的 Majorana（2 个 = 1 个 Dirac 费米子），非 D 类单个 Majorana。")
        print("  → 这是框架独有公式候选：标准 spinless（T²=+1）到不了 DIII，只有涌现 T'²=−1 能。")
    elif classes == {"D（单 Majorana）"}:
        print("  配对破 T'，BdG 回到 D 类（标准 spinless p+ip），接口非独有。")
    else:
        print(f"  混合/未定：{classes}，需进一步分析。")

    summary = {
        "question": "does the intrinsic Kramers T'=JK (T'^2=-1) survive inter-sublattice pairing, making the BdG class DIII instead of D?",
        "criterion": "T' survives iff J anti-commutes with Delta (J Delta J^-1 = -Delta, since Delta*=-Delta)",
        "results": results,
        "conclusion": (
            "if J anti-commutes with Delta, BdG is DIII (Kramers-pair Majorana at vortex); "
            "if J commutes, pairing breaks T' and BdG is D (single Majorana, standard spinless p+ip)"
        ),
    }
    out = ROOT / "experiments" / "exp_kramers_pairing_class_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

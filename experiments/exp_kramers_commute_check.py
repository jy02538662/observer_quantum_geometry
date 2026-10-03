"""决定性验证：涌现 J（谷 Kramers）与「spinless 子格间配对」的对易/反对易关系。

核心问题（探索独有公式的关键分岔）：
  涌现 T'=JK 在 BdG 存活 ⟺ J Δ Jᵀ = −Δ（Δ 纯虚，J 实反对称正交 ⟹ J⁻¹=Jᵀ=−J）。
  反对易（JΔJᵀ=−Δ）⟹ DIII（Kramers 对 Majorana）；对易 ⟹ D（单 Majorana）。

关键物理（泡利统计 + 子格结构）：
  spinless 费米子 + 2 子格 ⟹ 费米统计 ⟹ 配对 Δ=−Δᵀ。
  反对称性可来自「子格间」（iσ_y 子格反对称）而非「动量奇」。
  若配对是「动量偶（s 波）」→ 谷对称 → 应保 T'（DIII）；
  若配对是「动量奇（p 波）」→ 谷奇 → 破 T'（D）。

本实验直接算 J Δ Jᵀ 对几种配对（原胞内 / 键 s 波 / 键 p 波），看谁保 T'。

Code: `py -m experiments.exp_kramers_commute_check`
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
            return None
        sub = evec[:, i:j]
        for k in range(0, m, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J


def pairing_intracell(L, Delta0=0.5):
    """原胞内 A-B 配对（沿 x，动量无关 s 波）。"""
    N = L * L
    D = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for y in range(L):
        for x in range(0, L, 2):
            i = idx(x, y)
            j = idx(x + 1, y)
            D[i, j] += 1j * Delta0
            D[j, i] -= 1j * Delta0
    return D


def pairing_bond_s(L, Delta0=0.5):
    """键 s 波：所有 A-B 近邻（x,y 方向）同号 iΔ₀。"""
    N = L * L
    D = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            for (dx, dy) in [(1, 0), (0, 1)]:
                j = idx(x + dx, y + dy)
                D[i, j] += 1j * Delta0 / 2.0
                D[j, i] -= 1j * Delta0 / 2.0
    return D


def pairing_bond_p(L, Delta0=0.5):
    """键 p 波（手征 p+ip）：x 键 Δ₀、y 键 iΔ₀（动量奇）。"""
    N = L * L
    D = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            D[i, j] += Delta0 / 2.0
            D[j, i] -= Delta0 / 2.0
            j = idx(x, y + 1)
            D[i, j] += 1j * Delta0 / 2.0
            D[j, i] -= 1j * Delta0 / 2.0
    return D


def check_commute(J, Delta, name):
    JdJ = J @ Delta @ J.T  # J⁻¹ = Jᵀ = −J（实反对称正交）
    anti = np.allclose(JdJ, -Delta, atol=1e-8)
    comm = np.allclose(JdJ, +Delta, atol=1e-8)
    # 也看「混合度」：JdJ 在 Δ 方向的分量
    norm = np.linalg.norm(Delta)
    overlap_anti = np.real(np.vdot(-Delta, JdJ)) / (norm * np.linalg.norm(JdJ)) if norm > 0 else 0
    overlap_comm = np.real(np.vdot(+Delta, JdJ)) / (norm * np.linalg.norm(JdJ)) if norm > 0 else 0
    cls = "DIII（反对易，T'存活）" if anti else ("D（对易，T'破）" if comm else "混合")
    print(f"  {name}:")
    print(f"    JΔJᵀ=−Δ（反对易/DIII）: {anti},  JΔJᵀ=+Δ（对易/D）: {comm}")
    print(f"    与 −Δ 重叠 = {overlap_anti:+.3f}, 与 +Δ 重叠 = {overlap_comm:+.3f} → {cls}")
    return {"anti": bool(anti), "comm": bool(comm), "overlap_anti": round(float(overlap_anti), 3), "overlap_comm": round(float(overlap_comm), 3), "class": cls}


def main():
    print("=== 涌现 J 与各类 spinless 配对的对易关系 ===")
    print()

    L = 8
    H = pi_flux(L)
    J = find_J(H)
    print(f"[L={L}] J 性质（J²=−1, JD=DJ）: {np.allclose(J@J, -np.eye(L*L), atol=1e-9) and np.allclose(J@H, H@J, atol=1e-9)}")
    print()

    results = {}
    for name, D in [("原胞内（s 波谷对称）", pairing_intracell(L)),
                    ("键 s 波", pairing_bond_s(L)),
                    ("键 p+ip（动量奇）", pairing_bond_p(L))]:
        results[name] = check_commute(J, D, name)
        print()

    print("=== 结论 ===")
    print("  关键：框架的「子格间 iσ_y 配对」是 s 波（动量偶，谷对称）还是 p 波（动量奇，谷奇）？")
    print("  若 s 波谷对称 → 与 J 反对易 → DIII 类（Kramers 对 Majorana）= 独有公式。")
    print("  若 p 波谷奇 → 与 J 对易 → D 类（单 Majorana）= 接口非独有。")

    summary = {"L": L, "results": results}
    out = ROOT / "experiments" / "exp_kramers_commute_check_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

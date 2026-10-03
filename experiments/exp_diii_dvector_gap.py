"""攻分岔第七步（最终）：谷三重态 d 矢量配对在实空间完整模型的完全 gap 验证。

符号已坐实（exp_diii_dvector_final）：
  Δ = σ_y ⊗ (d·τ)，费米统计反对称 ✓，涌现 T' 存活（DIII）✓。

现在数值验证完全 gap。关键是 d(k) 动量奇覆盖 Dirac 点所有方向。

实空间构造：π 磁通完整模型（N=L²），涌现 J（find_J），谷三重态 τ_x/τ_z（已精确构造）。
配对 Δ = σ_y(子格间) ⊗ (d·τ)(谷三重态 × 动量奇)。

但实空间里 σ_y 和 τ 已折叠进 N 维。正确的实空间配对 = 用 τ_x/τ_z（非局域谷算符）
构造 Δ = iτ_x·P + iτ_z·P 型？前面发现 τ·P 矩阵乘不对易。

关键修正：τ 和 P 是不同自由度（谷 vs 轨道），应以「反对称化张量积」组合。
正确构造：配对 Δ_ij 在 BdG 里，需同时满足：
  (1) 费米统计 Δ=−Δᵀ；
  (2) 时间反演偶（T' 存活）；
  (3) 动量奇（覆盖 Dirac 点）。

让我用「对易子方案」：Δ = (τ_z P + P τ_z)/2 或 τ 与 P 的反对称化/对称化组合，
找一个既反对称（费米统计）又 T' 存活（DIII）又完全 gap 的配对。

Code: `py -m experiments.exp_diii_dvector_gap`
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


def tau_xz(H, tol=1e-8):
    ev, evec = np.linalg.eigh(H)
    N = len(ev)
    tx = np.zeros((N, N))
    tz = np.zeros((N, N))
    i = 0
    while i < N:
        j = i
        while j < N and abs(ev[j] - ev[i]) < tol:
            j += 1
        m = j - i
        sub = evec[:, i:j]
        for k in range(0, m, 2):
            tx += np.outer(sub[:, k + 1], sub[:, k]) + np.outer(sub[:, k], sub[:, k + 1])
            tz += np.outer(sub[:, k], sub[:, k]) - np.outer(sub[:, k + 1], sub[:, k + 1])
        i = j
    return tx, tz


def p_ip_pairing(L, D0):
    N = L * L
    P = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            P[i, j] += D0 / 2.0
            P[j, i] -= D0 / 2.0
            j = idx(x, y + 1)
            P[i, j] += 1j * D0 / 2.0
            P[j, i] -= 1j * D0 / 2.0
    return P


def build_and_measure(L, D0, pairing_kind):
    N = L * L
    H = pi_flux(L)
    J = find_J(H)
    tx, tz = tau_xz(H)
    P = p_ip_pairing(L, D0)

    if pairing_kind == "tauz_anticomm":
        # Δ = (τ_z P + P τ_z)/2（对称化，让 Δ 与 τ_z 对易，从而与 J 反对易 → T' 存活）
        D = (tz @ P + P @ tz) / 2.0
    elif pairing_kind == "tauz_comm":
        D = (tz @ P - P @ tz) / 2.0
    elif pairing_kind == "pure_p":
        D = P

    # 反对称化（费米统计）
    D = (D - D.T) / 2.0

    Z = np.zeros((N, N), dtype=complex)
    Hb = np.block([[H, D], [D.conj().T, -H]]).astype(complex)
    ev = np.linalg.eigvalsh(Hb)
    gap = np.min(np.abs(ev))

    # T' 存活
    Jbdg = np.block([[J, np.zeros((N, N))], [np.zeros((N, N)), J]])
    Tprime = Jbdg @ Hb.conjugate() @ Jbdg.T
    rel = np.linalg.norm(Tprime - Hb) / np.linalg.norm(Hb)
    return gap, rel, np.linalg.norm(D + D.T)


def main():
    print("=== 谷三重态 d 矢量配对：完全 gap 验证 ===")
    print()

    L = 8
    for kind in ("pure_p", "tauz_anticomm", "tauz_comm"):
        gap, rel, asym = build_and_measure(L, 0.5, kind)
        print(f"  {kind}: min gap={gap:.5f}, T'存活(rel={rel:.3f}), Δ反对称残差={asym:.2e}")

    summary = {"L": L}
    out = ROOT / "experiments" / "exp_diii_dvector_gap_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

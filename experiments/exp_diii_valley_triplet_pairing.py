"""攻分岔第五步：构造「谷三重态 d 矢量」配对，验证 T' 存活（DIII）+ 完全 gap。

核心（3He-B 机制框架对应）：
  涌现 J 对谷三重态 τ_x/τ_z 反对易 ⟹ τ 在 T' 下翻号。
  动量奇 d(k) 在时间反演下也翻号。两者乘积 T' 偶 ⟹ DIII。
  且动量奇配对覆盖所有 Dirac 点方向 ⟹ 完全 gap。

配对构造（实空间 N×N，spinless 费米统计 Δ=−Δᵀ）：
  标准 p+ip（完全 gap，D 类）：x 键 Δ₀、y 键 iΔ₀（子格间，动量奇，反对称）。
  谷三重态版本：Δ = τ ⊗ P，其中 τ 谷三重态（对称）、P 动量奇键配对（反对称）。
  但矩阵乘 τP 要反对称 ⟹ 需 [τ, P]=0。

先检查 τ_x/τ_z 与 p+ip 配对矩阵 P 的对易关系，确定正确构造。

Code: `py -m experiments.exp_diii_valley_triplet_pairing`
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
    """标准 p+ip 配对（反对称，动量奇，完全 gap，D 类）。"""
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


def main():
    print("=== 谷三重态 d 矢量配对：构造 + T' 存活 + 完全 gap ===")
    print()

    L = 8
    H = pi_flux(L)
    N = L * L
    J = find_J(H)
    tx, tz = tau_xz(H)

    P = p_ip_pairing(L, 0.5)

    print("1. τ_x/τ_z 与 p+ip 配对矩阵 P 的对易关系")
    comm_tx = np.allclose(tx @ P, P @ tx, atol=1e-8)
    comm_tz = np.allclose(tz @ P, P @ tz, atol=1e-8)
    anti_tx = np.allclose(tx @ P, -P @ tx, atol=1e-8)
    anti_tz = np.allclose(tz @ P, -P @ tz, atol=1e-8)
    print(f"   [τ_x, P]=0: {comm_tx}, {{τ_x, P}}=0: {anti_tx}")
    print(f"   [τ_z, P]=0: {comm_tz}, {{τ_z, P}}=0: {anti_tz}")
    print()

    # 构造谷三重态配对：Δ = τ ⊗ P 的几种组合，找反对称且保 T' 的
    print("2. 构造候选配对（需费米统计反对称 Δ=−Δᵀ）")
    candidates = {}
    # τ_z 与 P 反对易 ⟹ Δ = τ_z @ P 反对称？检查
    for name, T, tag in [("τ_z·P", tz, "τz"), ("τ_x·P", tx, "tx")]:
        D = T @ P
        anti = np.allclose(D, -D.T, atol=1e-8)
        candidates[name] = D
        print(f"   Δ = {name}: 反对称(Δ=−Δᵀ) = {anti}")
    print()

    # 3. 验证 T' 存活 + gap
    print("3. T' 存活 + 完全 gap 验证")
    Z = np.zeros((N, N), dtype=complex)
    Jbdg = np.block([[J, np.zeros((N, N))], [np.zeros((N, N)), J]])
    for name, D in candidates.items():
        Hb = np.block([[H, D], [D.conj().T, -H]]).astype(complex)
        # T' 存活
        Tprime = Jbdg @ Hb.conjugate() @ Jbdg.T
        rel = np.linalg.norm(Tprime - Hb) / np.linalg.norm(Hb)
        # gap
        ev = np.linalg.eigvalsh(Hb)
        gap = np.min(np.abs(ev))
        print(f"   Δ={name}: T'存活(rel={rel:.2e}) → {'DIII' if rel<1e-8 else 'D'}, min gap={gap:.5f}")

    summary = {"L": L, "comm_tx_P": bool(comm_tx), "comm_tz_P": bool(comm_tz), "anti_tx_P": bool(anti_tx), "anti_tz_P": bool(anti_tz)}
    out = ROOT / "experiments" / "exp_diii_valley_triplet_pairing_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

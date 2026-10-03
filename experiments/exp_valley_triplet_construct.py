"""攻分岔第四步：从涌现 J 精确构造谷三重态 τ 算符。

之前用 (−1)^{x+y} 猜测「谷间交错」是错的（min gap 仍 0）。正确做法：
从涌现 J（Kramers，J²=−1, JD=DJ, 实反对称正交）精确构造「谷三重态」。

谷三重态 τ_x, τ_y, τ_z 的定义（在谷空间 = Kramers 对空间）：
  - 厄米 τ_i²=I；
  - [τ_i, D]=0（内部对称，与哈密顿量对易）；
  - {τ_i, J}=0（与 J 反对易 ⟹ 在涌现 T'=JK 下翻号，这是 3He-B 机制的框架对应）。

其中 τ_y = J 本身（实反对称）。τ_x, τ_z 是实对称、与 J 反对易、与 D 对易。

构造方法：J 在每个 Kramers 对本征子空间 {ψ, Jψ} 是 [[0,1],[-1,0]]。
τ_x = 把 ψ↔Jψ 的实对称（= 每个子空间的 [[0,1],[1,0]] 或 [[1,0],[0,-1]] 组合）。
先数值构造，再验证四条性质。

Code: `py -m experiments.exp_valley_triplet_construct`
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
    """涌现 J：实反对称正交，J²=−1，JD=DJ。"""
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


def construct_tau_xz(H, J, tol=1e-8):
    """构造谷三重态的 τ_x, τ_z：实对称、τ²=I、[τ,D]=0、{τ,J}=0。

    在每个 Kramers 对本征子空间 {ψ, Jψ} 里，J=[[0,1],[-1,0]]。
    τ_x = [[0,1],[1,0]]（实对称，反对易 J），τ_z = [[1,0],[0,-1]]（实对称，反对易 J）。
    """
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
        sub = evec[:, i:j]  # N x m，本征子空间（m 偶，Kramers 对）
        for k in range(0, m, 2):
            # 在 {sub[:,k], sub[:,k+1]} 这 2 维里放 τ_x=[[0,1],[1,0]], τ_z=[[1,0],[0,-1]]
            # 但需确保 J 在这 2 维是 [[0,1],[-1,0]]（即 sub[:,k+1] = J sub[:,k]）
            tx += np.outer(sub[:, k + 1], sub[:, k]) + np.outer(sub[:, k], sub[:, k + 1])
            tz += np.outer(sub[:, k], sub[:, k]) - np.outer(sub[:, k + 1], sub[:, k + 1])
        i = j
    return tx, tz


def main():
    print("=== 从涌现 J 精确构造谷三重态 τ_x, τ_z ===")
    print()

    L = 8
    H = pi_flux(L)
    N = L * L
    J = find_J(H)
    tx, tz = construct_tau_xz(H, J)

    print("1. J 的性质（复核）")
    print(f"   J²=−1: {np.allclose(J@J, -np.eye(N), atol=1e-8)}")
    print(f"   JD=DJ: {np.allclose(J@H, H@J, atol=1e-8)}")
    print(f"   J 反对称: {np.allclose(J, -J.T, atol=1e-8)}")
    print()

    print("2. τ_x, τ_z 的四条性质")
    for name, T in [("τ_x", tx), ("τ_z", tz)]:
        herm = np.allclose(T, T.T, atol=1e-8)
        sq = np.allclose(T @ T, np.eye(N), atol=1e-8)
        commD = np.allclose(T @ H, H @ T, atol=1e-8)
        antiJ = np.allclose(T @ J + J @ T, np.zeros((N, N)), atol=1e-8)
        print(f"   {name}: 厄米={herm}, τ²=I={sq}, [τ,D]=0={commD}, {{τ,J}}=0={antiJ}")
    print()

    print("3. τ_x, τ_z 对易关系（应 {τ_x, τ_z}=0，构成谷 SU(2)）")
    print(f"   {{τ_x, τ_z}}=0: {np.allclose(tx@tz + tz@tx, np.zeros((N,N)), atol=1e-8)}")
    print(f"   τ_x τ_z = iτ_y(=iJ)? 检查 τ_x τ_z 与 J 关系")

    summary = {
        "J_ok": bool(np.allclose(J @ J, -np.eye(N), atol=1e-8) and np.allclose(J @ H, H @ J, atol=1e-8)),
        "tau_x": {"hermitian": bool(np.allclose(tx, tx.T, atol=1e-8)), "square_I": bool(np.allclose(tx@tx, np.eye(N), atol=1e-8)),
                  "comm_D": bool(np.allclose(tx@H, H@tx, atol=1e-8)), "anti_J": bool(np.allclose(tx@J + J@tx, np.zeros((N,N)), atol=1e-8))},
        "tau_z": {"hermitian": bool(np.allclose(tz, tz.T, atol=1e-8)), "square_I": bool(np.allclose(tz@tz, np.eye(N), atol=1e-8)),
                  "comm_D": bool(np.allclose(tz@H, H@tz, atol=1e-8)), "anti_J": bool(np.allclose(tz@J + J@tz, np.zeros((N,N)), atol=1e-8))},
    }
    out = ROOT / "experiments" / "exp_valley_triplet_construct_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

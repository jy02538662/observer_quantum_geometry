"""
下一刀：构造「局域实反对称」DIII 配对（谷奇/d 波）。

上一刀：Δ=J（实，非局域 Kramers 算符）→ DIII（2e-15）；Δ=iJ（虚）→ D；s 波/p+ip（虚）→ D。
判据：DIII 配对是「实」（T'-偶，JΔ*Jᵀ=+Δ）。

这一刀：测多个「局域实反对称」候选配对，找 DIII 相容的那个。

候选（都实、都反对称 Δ=−Δᵀ）：
  ① 实 p 波（x 键 +Δ₀、y 键 −Δ₀）
  ② 实 p 波带交错（x 键 +Δ₀(−1)^y、y 键 −Δ₀(−1)^x）——谷奇候选
  ③ 交错实配对（x 键 +Δ₀(−1)^{x+y}）

Code: `py -m experiments.verify_DIII_local_pairing`
"""
import numpy as np
from experiments._common import report


def build_H(N, m):
    n = N * N
    H = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = (x % N) * N + ((y + 1) % N)
            H[i, j] -= (-1.0) ** x
            H[j, i] -= (-1.0) ** x
            H[i, i] += m * (-1.0) ** ((x + y) % 2)
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
        m_ = j - i
        if m_ % 2 != 0:
            return None, f"奇重数 {m_}"
        sub = evec[:, i:j]
        for k in range(0, m_, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J, None


def check_DIII(H, J, Delta, n):
    Z = np.zeros((n, n), dtype=complex)
    H_bdg = np.block([[H, Delta], [Delta.conj().T, -H]])
    Jbdg = np.block([[J, np.zeros((n, n))], [np.zeros((n, n)), J]])
    T_prime = Jbdg @ H_bdg.conjugate() @ Jbdg.T
    return np.linalg.norm(T_prime - H_bdg) / np.linalg.norm(H_bdg)


def real_p_wave(N, Delta0=0.5):
    """实 p 波：x 键 +Δ₀、y 键 −Δ₀。"""
    n = N * N
    D = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            D[i, j] += Delta0
            D[j, i] -= Delta0
            j = (x % N) * N + ((y + 1) % N)
            D[i, j] -= Delta0
            D[j, i] += Delta0
    return D


def real_p_wave_staggered(N, Delta0=0.5):
    """实 p 波带交错（谷奇候选）：x 键 +Δ₀(−1)^y、y 键 −Δ₀(−1)^x。"""
    n = N * N
    D = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            ph = (-1.0) ** y
            D[i, j] += Delta0 * ph
            D[j, i] -= Delta0 * ph
            j = (x % N) * N + ((y + 1) % N)
            ph = (-1.0) ** x
            D[i, j] -= Delta0 * ph
            D[j, i] += Delta0 * ph
    return D


def staggered_real(N, Delta0=0.5):
    """交错实配对：x 键 +Δ₀(−1)^{x+y}。"""
    n = N * N
    D = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            ph = (-1.0) ** ((x + y) % 2)
            D[i, j] += Delta0 * ph
            D[j, i] -= Delta0 * ph
    return D


def run():
    R = {}
    N = 8
    n = N * N
    m = 0.1
    H = build_H(N, m)
    J, err = find_J(H)
    if J is None:
        R["note"] = err
        return R

    candidates = {
        "real_p_wave": real_p_wave(N),
        "real_p_wave_staggered": real_p_wave_staggered(N),
        "staggered_real": staggered_real(N),
    }
    for name, D in candidates.items():
        diff = check_DIII(H, J, D, n)
        R[name] = {
            "rel_diff": float(f"{diff:.2e}"),
            "DIII": diff < 1e-8,
            "note": f"{name} → {'DIII' if diff < 1e-8 else 'D'}",
        }

    R["verdict"] = {
        "goal": "找局域实反对称配对，满足 JΔ*Jᵀ=+Δ（DIII）",
        "next": "若都破 T'，则需构造显式谷（τ）结构：τ_z = 谷算符，配对 = 谷奇实反对称",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_DIII_local_pairing")

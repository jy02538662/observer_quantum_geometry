"""实空间严格验证：谷奇 d 波（AFM）vs 谷偶 s 波（FM）的 D vs DIII 分类。

核心（sympy 精确已证，这里在实空间完整模型钉死）：
  DIII 配对 = 谷奇（谷三重态 τz，连接 Kramers 对 k↔k+(π,π)）
  D 配对 = 谷偶（谷单态 τI，谷内配对）

关键翻转：之前「谷偶=DIII、谷奇=D」搞反了。正确是「谷偶=D、谷奇=DIII」。
物理原因：DIII 是 Kramers 对（谷对）配对，需谷间（谷奇）配对连接两个 Kramers 伙伴；
谷内（谷偶）配对不连接 Kramers 对，反而破 T'。

实空间构造：
  谷间配对 = 连接 k 和 k+(π,π) = 实空间乘交错因子 (−1)^{x+y}（因 e^{i(π,π)·r}=(−1)^{x+y}）
  d 波 = cos kx − cos ky = x 键 +Δ₀、y 键 −Δ₀

候选（都在交错质量 m 下）：
  ① Δ=J（谷 Kramers 算符，非局域谷间配对）——基准，应 DIII
  ② 交错 d 波（谷奇）：x 键 +Δ₀(−1)^{x+y}、y 键 −Δ₀(−1)^{x+y}
  ③ 非交错 d 波（谷偶）：x 键 +Δ₀、y 键 −Δ₀
  ④ 交错 s 波（谷奇，无 d 波角动量）：x/y 键 +Δ₀(−1)^{x+y}

Code: `py -m experiments.verify_DIII_valley_parity`
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
    H_bdg = np.block([[H, Delta], [Delta.conj().T, -H]])
    Jbdg = np.block([[J, np.zeros((n, n))], [np.zeros((n, n)), J]])
    T_prime = Jbdg @ H_bdg.conjugate() @ Jbdg.T
    return np.linalg.norm(T_prime - H_bdg) / np.linalg.norm(H_bdg)


def staggered_d_wave(N, Delta0=0.5):
    """谷奇 d 波：x 键 +Δ₀(−1)^{x+y}、y 键 −Δ₀(−1)^{x+y}（交错 + d 波角动量）。"""
    n = N * N
    D = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            ph = (-1.0) ** ((x + y) % 2)
            j = ((x + 1) % N) * N + (y % N)
            D[i, j] += Delta0 * ph
            D[j, i] -= Delta0 * ph
            j = (x % N) * N + ((y + 1) % N)
            D[i, j] -= Delta0 * ph
            D[j, i] += Delta0 * ph
    return D


def plain_d_wave(N, Delta0=0.5):
    """谷偶 d 波：x 键 +Δ₀、y 键 −Δ₀（无交错因子）。"""
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


def staggered_s_wave(N, Delta0=0.5):
    """谷奇 s 波（交错，无 d 波角动量）：x/y 键 +Δ₀(−1)^{x+y}。"""
    n = N * N
    D = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            ph = (-1.0) ** ((x + y) % 2)
            j = ((x + 1) % N) * N + (y % N)
            D[i, j] += Delta0 * ph
            D[j, i] -= Delta0 * ph
            j = (x % N) * N + ((y + 1) % N)
            D[i, j] += Delta0 * ph
            D[j, i] -= Delta0 * ph
    return D


def run():
    R = {}
    for N in [8, 12]:
        n = N * N
        m = 0.1
        H = build_H(N, m)
        J, err = find_J(H)
        if J is None:
            R[f"N{N}"] = {"note": err}
            continue

        res = {}
        # ① 基准 Δ=J（非局域谷间）
        res["Delta_J（谷间，非局域）"] = float(f"{check_DIII(H, J, J, n):.2e}")
        # ② 谷奇 d 波（交错 d 波）
        res["谷奇 d 波（交错 d）"] = float(f"{check_DIII(H, J, staggered_d_wave(N), n):.2e}")
        # ③ 谷偶 d 波（无交错）
        res["谷偶 d 波（无交错）"] = float(f"{check_DIII(H, J, plain_d_wave(N), n):.2e}")
        # ④ 谷奇 s 波（交错 s）
        res["谷奇 s 波（交错）"] = float(f"{check_DIII(H, J, staggered_s_wave(N), n):.2e}")

        R[f"N{N}"] = {
            "rel_diff（<1e-8 = DIII）": res,
            "note": "Δ=J 应 DIII（基准）；谷奇 d 波若 <1e-8 则 DIII 钉死",
        }

    R["verdict"] = {
        "criterion": "rel_diff < 1e-8 ⟹ T'-偶 ⟹ DIII（保 Kramers）；否则 D（破 T'）",
        "expectation": "谷奇（交错，连接 Kramers 对）→ DIII；谷偶（无交错，谷内）→ D",
        "physics": "DIII 是 Kramers 对配对，需谷间（谷奇）连接；谷内（谷偶）不连接 Kramers 对反而破 T'",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_DIII_valley_parity")

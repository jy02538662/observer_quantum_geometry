"""
找 DIII 相容的配对对称性：测多个候选配对，看哪个满足 [T', H_BdG]=0。

判据（AZ 分类）：T' Δ T'^{-1} = Δ（DIII，Kramers 对 Majorana）vs = −Δ（D，单 Majorana）。

候选配对：
  ① s 波（原胞内子格间 iσ_y）：Δ_ij = iΔ₀（磁单位胞 (2j,y)-(2j+1,y)）→ 之前测出破 T'（D）
  ② 手征 p+ip（键方向）：x 键 Δ₀、y 键 iΔ₀（动量依赖 sin k_x + i sin k_y）
  ③ Δ = J（Kramers 算符本身当配对，toy 验证判据）——若 J J* Jᵀ = J 则 DIII 相容

Code: `py -m experiments.verify_DIII_pairing`
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


def s_wave_pairing(N, Delta0=0.5):
    n = N * N
    D = np.zeros((n, n), dtype=complex)
    for y in range(N):
        for x in range(0, N, 2):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            D[i, j] += 1j * Delta0
            D[j, i] -= 1j * Delta0
    return D


def pip_pairing(N, Delta0=0.5):
    """手征 p+ip：x 键 Δ₀、y 键 iΔ₀。"""
    n = N * N
    D = np.zeros((n, n), dtype=complex)
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            D[i, j] += Delta0
            D[j, i] -= Delta0
            j = (x % N) * N + ((y + 1) % N)
            D[i, j] += 1j * Delta0
            D[j, i] -= 1j * Delta0
    return D


def check_DIII(H, J, Delta, n):
    Z = np.zeros((n, n), dtype=complex)
    H_bdg = np.block([[H, Delta], [Delta.conj().T, -H]])
    Jbdg = np.block([[J, np.zeros((n, n))], [np.zeros((n, n)), J]])
    T_prime = Jbdg @ H_bdg.conjugate() @ Jbdg.T
    return np.linalg.norm(T_prime - H_bdg) / np.linalg.norm(H_bdg)


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

    # ① s 波
    D_s = s_wave_pairing(N)
    R["s_wave"] = {
        "rel_diff": float(f"{check_DIII(H, J, D_s, n):.2e}"),
        "DIII": check_DIII(H, J, D_s, n) < 1e-8,
        "note": "s 波（iσ_y 原胞内）→ 若 rel_diff 大则破 T'（D）",
    }

    # ② 手征 p+ip
    D_pip = pip_pairing(N)
    R["p_ip"] = {
        "rel_diff": float(f"{check_DIII(H, J, D_pip, n):.2e}"),
        "DIII": check_DIII(H, J, D_pip, n) < 1e-8,
        "note": "手征 p+ip（x 键 Δ₀、y 键 iΔ₀）",
    }

    # ③ Δ = J（实，toy 验证判据）——Δ=J 是 T'-偶（DIII），Δ=iJ 是 T'-奇（D）
    for label, Dj in [("Delta_J_real", J), ("Delta_iJ_imag", 1j * J)]:
        kind = "J（实，T'偶）" if np.allclose(Dj, J) else "iJ（虚，T'奇）"
        R[label] = {
            "rel_diff": float(f"{check_DIII(H, J, Dj, n):.2e}"),
            "DIII": check_DIII(H, J, Dj, n) < 1e-8,
            "note": f"Δ={kind}——验证 T'ΔT'⁻¹ 判据",
        }

    R["verdict"] = {
        "criterion": "T' Δ T'^{-1} = Δ（DIII）vs −Δ（D）",
        "next": "若 s 波破 T'、p+ip 也破、Δ=J 破/不破，则需构造「谷奇单态 Δ=iΔ₀ σ_y τ_z」或「d 波」",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_DIII_pairing")

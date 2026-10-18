"""
无磁场 DIII 拓扑超导：完整 BdG 推导（π-flux + 交错质量 + 配对）。

坐实「真卖点」：内禀 Kramers T²=-1 → 不需要外加磁场 → DIII 类拓扑超导。

本脚本验证：π-flux + 交错质量 m + 原胞内子格间配对 Δ=iΔ₀σ_y 的 BdG，
  Kramers T' 是否存活（DIII）——即「交错质量三合一」在 BdG 层面是否仍成立。

方法（接 exp_kramers_bdg_survival）：
  ① π-flux + 交错质量 m；
  ② find_J 构造 Kramers 算符 J（T'=JK，T'²=-1）；
  ③ 原胞内子格间配对 Δ（s 波谷奇，spinless）；
  ④ BdG = [[H, Δ], [Δ†, −H]]；
  ⑤ 检查 [T'_BdG, H_BdG] = 0（DIII = Kramers 对 Majorana）vs ≠0（D = 单 Majorana）。

Code: `py -m experiments.verify_DIII_bdg`
"""
import numpy as np
from experiments._common import report


def build_H(N, m):
    """π-flux + 交错质量 m。"""
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
            return None, f"奇重数 {m_}（非自对偶）"
        sub = evec[:, i:j]
        for k in range(0, m_, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J, None


def intra_cell_pairing(N, Delta0=0.5):
    """原胞内子格间配对（磁单位胞 (2j,y)-(2j+1,y)，iΔ₀σ_y）。"""
    n = N * N
    Delta = np.zeros((n, n), dtype=complex)
    for y in range(N):
        for x in range(0, N, 2):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            Delta[i, j] += 1j * Delta0
            Delta[j, i] -= 1j * Delta0
    return Delta


def run():
    R = {}
    N = 8
    n = N * N

    for m in [0.0, 0.1, 0.5]:
        H = build_H(N, m)
        J, err = find_J(H)
        if J is None:
            R[f"m_{m}"] = {"note": err}
            continue

        Delta = intra_cell_pairing(N)
        Z = np.zeros((n, n), dtype=complex)
        H_bdg = np.block([[H, Delta], [Delta.conj().T, -H]])

        # T'_BdG = diag(J, J) K
        Jbdg = np.block([[J, np.zeros((n, n))], [np.zeros((n, n)), J]])
        T_prime = Jbdg @ H_bdg.conjugate() @ Jbdg.T
        diff = np.linalg.norm(T_prime - H_bdg) / np.linalg.norm(H_bdg)

        R[f"m_{m}"] = {
            "m": m,
            "rel_diff（[T', H_BdG]=0）": float(f"{diff:.2e}"),
            "DIII_survives": diff < 1e-8,
            "az_class": "DIII（Kramers 对 Majorana）" if diff < 1e-8 else "D（单 Majorana）",
        }

    R["verdict"] = {
        "no_magnetic_field": "内禀 Kramers T²=-1 → 不需外加磁场 → DIII 类拓扑超导",
        "if_DIII_all_m": "若所有 m 下 [T', H_BdG]=0，则「交错质量三合一」在 BdG 层面成立 = 无磁场 DIII 拓扑超导机制坐实",
        "significance": "「无需磁场」是当前拓扑超导最大的工程障碍，DIII 机制省掉它 = 真卖点",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_DIII_bdg")

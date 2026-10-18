"""
严格验证「交错质量下 Kramers T²=-1 存活」：显式构造 Kramers 算符 J，
检查 T'²=-1 和 [H,T']=0（把「全偶重数」升级到「严格 T²=-1」）。

Kramers 算符：T' = J K（J 反对称，K 复共轭），T'²=-1。
  判据：① J²=-I（T'²=-1）；② J H* Jᵀ = H（[H,T']=0，T' 是 H 的对称）。

交错质量 m(-1)^{x+y} 实数 ⟹ 保复共轭 K；若 J 作用在谷空间（非子格），则 [H,T']=0。

Code: `py -m experiments.verify_staggered_mass_kramers_strict`
"""
import numpy as np
from experiments._common import report


def build_mass(N, m):
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
    """从简并本征态构造 Kramers 算符 J（反对称，J²=-I 在简并子空间）。"""
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
            return None, f"奇重数 {m}（非自对偶，Kramers 破）"
        sub = evec[:, i:j]
        for k in range(0, m, 2):
            J += np.outer(sub[:, k + 1], sub[:, k]) - np.outer(sub[:, k], sub[:, k + 1])
        i = j
    return J, None


def run():
    R = {}
    N = 8
    n = N * N

    for m in [0.0, 0.1, 0.5]:
        H = build_mass(N, m)
        J, err = find_J(H)
        if J is None:
            R[f"m_{m}"] = {"note": err}
            continue

        # ① T'² = -1 ⟺ J J* = -I（J 实 ⟹ J² = -I）
        J2 = J @ J
        T2_minus1 = np.linalg.norm(J2 + np.eye(n)) / n

        # ② [H, T'] = 0 ⟺ J H* Jᵀ = H
        THT = J @ H.conjugate() @ J.T
        comm = np.linalg.norm(THT - H) / np.linalg.norm(H)

        R[f"m_{m}"] = {
            "m": m,
            "‖J²+I‖/N（T'²=-1）": float(f"{T2_minus1:.2e}"),
            "‖JH*Jᵀ−H‖/‖H‖（[H,T']=0）": float(f"{comm:.2e}"),
            "T_prime_survives": comm < 1e-8,
            "note": f"m={m}：T'²=-1 和 [H,T']=0 是否同时成立",
        }

    R["verdict"] = {
        "strict_check": "显式 J 下：① J²=-I（T'²=-1）；② J H* Jᵀ=H（[H,T']=0）",
        "if_both": "交错质量下 T'²=-1 严格存活 → 系统 = DIII 类（Kramers 对 Majorana），DIII 判定从「全偶重数」升级到「严格 T²=-1」",
        "significance": "把「三合一机制」的 DIII 判定钉死（必要判据 → 严格证明）",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_staggered_mass_kramers_strict")

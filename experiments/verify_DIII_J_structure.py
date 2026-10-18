"""
下一刀：从 J 的结构里挖出「局域谷奇」DIII 配对。

上一刀：Δ=J（全 Kramers 算符）→ DIII（2e-15）；naive 局域实配对全破。
这一刀：检查 J 的空间结构（哪些格点对 |i−j| 有非零元），
  找到 J 的「最短程结构」，据此构造局域 DIII 配对。

J = Σ_k (|ψ_{2k+1}⟩⟨ψ_{2k}| − |ψ_{2k}⟩⟨ψ_{2k+1}|)，连接 Kramers 对（谷对 k↔k+(π,π)）。

Code: `py -m experiments.verify_DIII_J_structure`
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

    # J 的空间结构：按 |dx|+|dy|（曼哈顿距离）分组，看哪些距离有非零元
    dist_strength = {}
    for i in range(n):
        xi, yi = i // N, i % N
        for j in range(n):
            xj, yj = j // N, j % N
            dx = min(abs(xi - xj), N - abs(xi - xj))
            dy = min(abs(yi - yj), N - abs(yi - yj))
            d = dx + dy
            dist_strength.setdefault(d, 0.0)
            dist_strength[d] += abs(J[i, j]) ** 2

    total = sum(dist_strength.values())
    dist_frac = {d: round(v / total, 4) for d, v in sorted(dist_strength.items())}
    R["J_distance_structure"] = {
        "manhattan_distance_fraction": dist_frac,
        "note": "J 的 |J_ij|² 按曼哈顿距离的分布——哪个距离主导",
    }

    # 取最短程主导距离，构造局域配对（只保留该距离的 J 分量）
    best_d = max(dist_frac, key=dist_frac.get)
    J_local = np.zeros((n, n))
    for i in range(n):
        xi, yi = i // N, i % N
        for j in range(n):
            xj, yj = j // N, j % N
            dx = min(abs(xi - xj), N - abs(xi - xj))
            dy = min(abs(yi - yj), N - abs(yi - yj))
            if dx + dy == best_d:
                J_local[i, j] = J[i, j]

    diff_local = check_DIII(H, J, J_local, n)
    R["truncated_J"] = {
        "kept_distance": best_d,
        "rel_diff": float(f"{diff_local:.2e}"),
        "DIII": diff_local < 1e-8,
        "note": f"只保留曼哈顿距离 {best_d} 的 J 分量 → 是否 DIII",
    }

    R["verdict"] = {
        "J_structure": "J 主导的曼哈顿距离 = 谷 Kramers 结构的空间尺度",
        "local_pairing": "若截断 J 给 DIII，则找到了局域谷奇配对；否则需动量空间显式谷算子",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_DIII_J_structure")

"""
查「交错质量下 Kramers T²=-1 是否存活」——本征值是否依然成对（全偶重数）。

Kramers 定理：时间反演 T'²=-1 ⟹ 每个本征值至少二重简并。
  交错质量 m(-1)^{x+y} 是实数（保复共轭 K），但破手征（子格 σ_z）。
  问题：它是否破 Kramers（T'²=-1）？判据 = 本征值是否依然全偶重数。

Code: `py -m experiments.verify_staggered_mass_kramers`
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


def multiplicity_check(evals, tol=1e-6):
    """检查本征值是否全偶重数（Kramers 成对）。"""
    evals = np.sort(evals)
    # 按容差分组
    groups = []
    current = [evals[0]]
    for e in evals[1:]:
        if abs(e - current[0]) < tol:
            current.append(e)
        else:
            groups.append(current)
            current = [e]
    groups.append(current)
    mults = [len(g) for g in groups]
    all_even = all(m % 2 == 0 for m in mults)
    return mults, all_even


def run():
    R = {}
    N = 8
    n = N * N

    for m in [0.0, 0.1, 0.5]:
        H = build_mass(N, m)
        evals = np.linalg.eigvalsh(H)
        mults, all_even = multiplicity_check(evals)
        R[f"m_{m}"] = {
            "m": m,
            "multiplicities": mults,
            "all_even（Kramers 存活）": all_even,
            "note": f"m={m} 时本征值重数是否全偶（Kramers T²=-1 是否存活）",
        }

    # 更细扫描：找 Kramers 破裂的临界 m
    fine_ms = [0.0, 0.05, 0.1, 0.2, 0.3]
    surv = {}
    for m in fine_ms:
        H = build_mass(N, m)
        evals = np.linalg.eigvalsh(H)
        _, all_even = multiplicity_check(evals)
        surv[m] = all_even
    R["fine_scan"] = {
        "m": fine_ms,
        "kramers_survives": surv,
        "note": "细扫描找 Kramers 是否随 m 破裂",
    }

    R["verdict"] = {
        "kramers_under_mass": "交错质量（实数，破手征）下 Kramers T²=-1 是否存活 = 本征值是否全偶重数",
        "if_survives": "若存活：交错质量 = 「破手征但保 Kramers」的开隙项 → 抑制 CDW + 保留吸引 + 保留 DIII 资产",
        "if_breaks": "若破裂：交错质量牺牲 Kramers，需换「保 T 的开隙项」（Kekulé/谷间，非 Haldane——Haldane 破 T）",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_staggered_mass_kramers")

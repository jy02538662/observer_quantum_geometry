"""
交错质量对照实验：验证「开单粒子隙 → 抑制 Peierls（|K| 减小）」这个普适机制。

物理：Peierls（K<0）在半填充（Dirac 点）最强，因为态密度为零、嵌套完美。
  加交错质量 m(-1)^{x+y} 开隙 Δ=2m，态密度不再为零 → Peierls 应被抑制（|K| 减小）。

但注意（用户诊断）：交错质量破坏手征（子格）对称性 = Kramers T²=-1 和 U<0（吸引）的根源，
  所以它是「对照实验」（验证机制），不是最终机制（可能同时削弱吸引）。

本脚本验证：
  A. 交错质量开隙（Δ = 2m）；
  B. |K| 随 m 增大而减小（Peierls 被抑制）；
  C. K 的符号是否保持负（吸引是否保留，还是被手征破缺削弱）。

Code: `py -m experiments.verify_staggered_mass_cdw`
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
            # 交错质量 m(-1)^{x+y}
            H[i, i] += m * (-1.0) ** ((x + y) % 2)
    return H


def dimer_dD(N):
    n = N * N
    dD = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            phase = (-1.0) ** ((x + y) % 2)
            dD[i, j] += phase
            dD[j, i] += phase
    dD /= np.linalg.norm(dD, "fro")
    return dD


def gap(H):
    return float(np.min(np.abs(np.linalg.eigvalsh(H))))


def stiffness(H, dD, nocc, eps=1e-3):
    E0 = np.sum(np.linalg.eigvalsh(H)[:nocc])
    Ep = np.sum(np.linalg.eigvalsh(H + eps * dD)[:nocc])
    Em = np.sum(np.linalg.eigvalsh(H - eps * dD)[:nocc])
    return (Ep - 2 * E0 + Em) / (eps**2)


def run():
    R = {}
    N = 8
    n = N * N
    nocc = n // 2
    dD = dimer_dD(N)

    m_list = [0.0, 0.1, 0.3, 0.5, 0.8, 1.0]
    gaps = []
    Ks = []
    for m in m_list:
        H = build_mass(N, m)
        gaps.append(round(gap(H), 4))
        Ks.append(round(stiffness(H, dD, nocc), 4))

    R["1_gap_and_K_vs_m"] = {
        "m（交错质量）": m_list,
        "Δ（能隙）": gaps,
        "K（二聚化刚度）": Ks,
        "gap_opens": gaps[-1] > gaps[0],
        "|K|_decreases": abs(Ks[-1]) < abs(Ks[0]),
        "K_stays_negative": all(k < 0 for k in Ks),
        "note": "开隙 → |K| 减小（Peierls 抑制）？K 保持负（吸引保留）？",
    }

    R["verdict"] = {
        "mechanism": "若「开隙 → |K| 减小」成立，则「开单粒子隙抑制 Peierls」是普适机制",
        "chiral_warning": "但交错质量破手征 → 可能削弱 U<0（吸引）——这是对照实验，非最终机制",
        "next": "若机制成立，下一步攻「正确 Haldane 项（开隙且保手征）」",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_staggered_mass_cdw")

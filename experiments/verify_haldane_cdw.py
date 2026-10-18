"""
算「Haldane 阻挫能否抑制 CDW（二聚化）而保持吸引（K<0）」。

物理：Peierls（K<0）= 吸引（负 U），想二聚化开隙变绝缘体。
  加 Haldane 型次近邻（复相位）＝ 阻挫，可能：
  ① 保持 K < 0（吸引保留）；
  ② 减小二聚化（CDW）能隙 Δ_CDW，让超导相竞争胜出。

本脚本验证：
  A. Haldane 次近邻开拓扑隙 Δ_Haldane（gap 随 t2 增大）；
  B. 二聚化刚度 K 随 t2 怎么变——是否保持负（吸引）但 |K| 减小（Peierls 被抑制）。

Code: `py -m experiments.verify_haldane_cdw`
"""
import numpy as np
from experiments._common import report


def build_haldane(N, t2, phi=np.pi / 2):
    n = N * N
    H = np.zeros((n, n), dtype=complex)

    def idx(x, y):
        return (x % N) * N + (y % N)

    # NN：π-flux
    for x in range(N):
        for y in range(N):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            H[i, j] -= (-1.0) ** x
            H[j, i] -= (-1.0) ** x

    # NNN：Haldane 型对角（复相位，破时间反演）
    for x in range(N):
        for y in range(N):
            i = idx(x, y)
            # (x+1, y+1) 与 (x-1, y-1)：+φ
            j = idx(x + 1, y + 1)
            H[i, j] += t2 * np.exp(1j * phi)
            H[j, i] += t2 * np.exp(-1j * phi)
            # (x+1, y-1) 与 (x-1, y+1)：−φ
            j = idx(x + 1, y - 1)
            H[i, j] += t2 * np.exp(-1j * phi)
            H[j, i] += t2 * np.exp(1j * phi)

    return H


def gap(H):
    """费米能级（0）处的能隙 = min |eigenvalue|。"""
    evals = np.linalg.eigvalsh(H)
    return float(np.min(np.abs(evals)))


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


def stiffness(H, dD, nocc, eps=1e-3):
    E0 = np.sum(np.linalg.eigvalsh(H)[:nocc]).real
    Ep = np.sum(np.linalg.eigvalsh(H + eps * dD)[:nocc]).real
    Em = np.sum(np.linalg.eigvalsh(H - eps * dD)[:nocc]).real
    return (Ep - 2 * E0 + Em) / (eps**2)


def run():
    R = {}
    N = 8
    n = N * N
    nocc = n // 2
    dD = dimer_dD(N)

    t2_list = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    gaps = []
    Ks = []
    for t2 in t2_list:
        H = build_haldane(N, t2)
        gaps.append(round(gap(H), 6))
        Ks.append(round(stiffness(H, dD, nocc), 4))

    R["1_gap_vs_t2"] = {
        "t2（Haldane 强度）": t2_list,
        "Δ_Haldane（拓扑隙）": gaps,
        "note": "Haldane 次近邻开拓扑隙（随 t2 增大）",
    }

    R["2_K_vs_t2"] = {
        "t2": t2_list,
        "K（二聚化刚度）": Ks,
        "K_stays_negative": all(k < 0 for k in Ks),
        "|K|_decreases": abs(Ks[-1]) < abs(Ks[0]),
        "note": "K 是否保持负（吸引）+ |K| 是否随 t2 减小（Peierls 被抑制）",
    }

    R["verdict"] = {
        "attractive_preserved": all(k < 0 for k in Ks),
        "cdw_suppressed": abs(Ks[-1]) < abs(Ks[0]),
        "conclusion": "若 K 保持负（吸引保留）且 |K| 减小（CDW 被抑制），则 Haldane 阻挫是「稳定 δD 而不杀超导」的正路",
        "next": "接着算超导 Tc vs CDW 能隙的竞争（哪个在低温赢）",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_haldane_cdw")

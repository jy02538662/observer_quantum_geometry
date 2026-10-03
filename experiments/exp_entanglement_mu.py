"""验证「掺杂 = 纠缠谱化学势」这个破墙猜测（第一步：纠缠谱连续 vs 离散）。

背景（承接 [[超导独有公式探索]] §九 + 讨论「μ 墙 = 系统/外界纠缠的产物」）：
  十八轮判负全失败，因为所有候选（断裂/绕数/拓扑荷）都是「离散」量，给不出「连续填充」。
  新猜测：μ = 子块 A 的纠缠谱化学势。对无能隙系统（π 磁通有 Dirac 点），
  约化密度矩阵 ρ_A 的本征值谱是连续的，天然带一个「费米能级/化学势」。

本实验：π 磁通 D（L×L 环面），spinless 费米子半满占据（负能态全占，Dirac 点正好在费米能级），
取「左半格」为子块 A，算约化密度矩阵 ρ_A = Tr_B|ψ><ψ| 的纠缠谱，
看它连续 vs 离散、有无自然化学势位置。

Code: `py -m experiments.exp_entanglement_mu`
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


def entanglement_spectrum(H, subA_idx):
    """基态（半满：负能态全占）的约化密度矩阵纠缠谱。"""
    ev, evec = np.linalg.eigh(H)
    N = len(ev)
    # 半满占据：负能态全占（费米能级在 0，Dirac 点）
    occ = ev < 0
    # 多体 Slater 行列式基态（占据态的 Slater 行列式）
    # 约化密度矩阵：对占据态 |ψ_k>，ρ_A = Σ_occ |ψ_k^A><ψ_k^A| 的关联矩阵
    # 对 free fermion，ρ_A 的本征值由关联矩阵 C_ij = <c_i† c_j>（i,j ∈ A）决定
    # C = Σ_occ ψ_k(i) ψ_k(j)*（i,j ∈ A），ρ_A 纠缠谱 = C 的本征值（占据概率）

    # 关联矩阵（子块 A 内）
    occ_vec = evec[:, occ]  # N × n_occ
    occ_A = occ_vec[subA_idx, :]  # n_A × n_occ
    C = occ_A @ occ_A.conj().T  # n_A × n_A 关联矩阵
    # 纠缠谱 = C 的本征值（0~1 之间的占据概率）
    spec = np.linalg.eigvalsh(C)
    return spec


def main():
    print("=== π 磁通子块纠缠谱：连续 vs 离散，有无化学势 ===")
    print()

    results = {}
    for L in (8, 12):
        H = pi_flux(L)
        N = L * L
        # 左半格为子块 A：x < L/2 的格点
        def idx(x, y):
            return (x % L) * L + (y % L)

        subA = [idx(x, y) for x in range(L // 2) for y in range(L)]
        spec = entanglement_spectrum(H, subA)
        spec = np.sort(spec)
        nA = len(subA)
        print(f"[L={L}, 子块 A 大小 n_A={nA}, 总 N={N}]")
        print(f"  纠缠谱（C 本征值，0~1 占据概率）范围: [{spec.min():.4f}, {spec.max():.4f}]")

        # 看分布：多少个本征值落在 (0,1) 内部（连续占据 vs 只有 0/1 离散）
        interior = np.sum((spec > 0.05) & (spec < 0.95))
        zero_one = np.sum(spec < 0.05) + np.sum(spec > 0.95)
        print(f"  内部(0.05~0.95)本征值数 = {interior}，接近0/1的数 = {zero_one}，总 = {nA}")
        print(f"  → {'连续分布（有中间占据）' if interior > 0 else '离散（只有0/1）'}")

        # 化学势位置：占据概率 = 0.5 处（若有）
        # 分位数分布
        qs = [0.1, 0.25, 0.5, 0.75, 0.9]
        qvals = [float(np.quantile(spec, q)) for q in qs]
        print(f"  分位数 {qs}: {[round(v,3) for v in qvals]}")
        print(f"  前 20 个本征值: {[round(v,4) for v in spec[:20]]}")
        print(f"  后 20 个本征值: {[round(v,4) for v in spec[-20:]]}")
        results[f"L{L}"] = {
            "nA": nA, "range": [float(spec.min()), float(spec.max())],
            "interior": int(interior), "zero_one": int(zero_one),
            "continuous": bool(interior > 0),
            "quantiles": {str(q): round(qvals[i], 3) for i, q in enumerate(qs)},
        }
        print()

    summary = {"question": "is the pi-flux entanglement spectrum continuous (with a natural chemical potential) or discrete?", "results": results}
    out = ROOT / "experiments" / "exp_entanglement_mu_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

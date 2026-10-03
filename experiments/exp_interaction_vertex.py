"""
四费米子相互作用顶点：Q[δD]（Hessian）重新解读为「δD 玻色场介导的相互作用」

背景（承接预印本 1.3 §4.2「涌现引力重新推导」的 Hessian 推导）：
  纯迹作用量 S = -αTr(D²) + βTr(D⁴) 在 π-flux D₀ 处的二阶变分
    Q[δD] = -αTr(δD²) + 4βTr(D₀²δD²) + 2βTr(D₀δD D₀δD)
  之前只把它当「Hessian 稳定性判据」（本征值 → 鞍点/极小）。
  本脚本重新解读：Q[δD] 是「键涨落 δD 作为玻色场」的有效二次作用量。

  物理对应（标准，非牵强）：
    δD_ij = 键 i-j 的涨落 = 玻色场（声子/磁振子）
    Q[δD] = 该玻色场的有效作用量（动能 + 质量项）
    本征值 ω_a = 各玻色子模式的质量²
    传播子 = 1/ω_a（软模 ω_a→0 给长程相互作用）
    δD 通过 Yukawa 耦合 S_int = g·Σ δD_ij c_i†c_j 耦合到费米子
    积分掉 δD ⟹ 四费米子 V_eff = g² · c†c · Q⁻¹ · c†c

  要回答的核心问题：
    软模（近零本征值）是「模长型（声子，键伸缩）」还是「相位型（磁振子，相位旋转）」？
    → 决定框架内生的是「配对相互作用（声子介导，吸引）」还是「密度相互作用（磁振子介导）」。

  基 {M_a} 的两部分（basis_matrices）：
    前 ne 个：对称实 M_ij = M_ji = 1（模长方向 δr_ij，键伸缩 = 声子）
    后 ne 个：反对称虚 M_ij = i, M_ji = -i（相位方向 δθ_ij，相位旋转 = 磁振子/规范场）

  ⚠️ 诚实边界（2026-10-03 运行后）：π-flux D₀ 是鞍点（Hessian 负本征值），
  「软模 = 玻色子」解读在鞍点处不成立（负质量² = 不稳定模）。本脚本的价值 =
  确认「四费米子解读需在真实基态（稳定极小）处做」，并读出鞍点下降方向的
  模长/相位混合结构。真实基态处（Hessian 正定）的正本征值才是物理玻色子。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


# ---- 内联的依赖函数（原 exp_gravity_mass_gap.py，已自包含，不 import 外部）----

def toroidal_D(n_per_dim, pi_flux=True):
    n = n_per_dim ** 2
    D = np.zeros((n, n), complex)
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            idx = n_per_dim * i + j
            jr = (j + 1) % n_per_dim
            D[idx, n_per_dim * i + jr] += 1.0
            D[n_per_dim * i + jr, idx] += 1.0
            ph = np.pi * j if pi_flux else 0.0
            D[idx, n_per_dim * ((i + 1) % n_per_dim) + j] += np.exp(1j * ph)
            D[n_per_dim * ((i + 1) % n_per_dim) + j, idx] += np.exp(-1j * ph)
    return D


def basis_matrices(n):
    ne = n * (n - 1) // 2
    Ms = []
    for i in range(n):
        for j in range(i + 1, n):
            M = np.zeros((n, n), complex); M[i, j] = 1.0; M[j, i] = 1.0; Ms.append(M)
    for i in range(n):
        for j in range(i + 1, n):
            M = np.zeros((n, n), complex); M[i, j] = 1j; M[j, i] = -1j; Ms.append(M)
    return np.asarray(Ms)


def simple_hessian(D0, alpha, beta):
    n = D0.shape[0]; D02 = D0 @ D0
    Ms = basis_matrices(n)
    LMs = (-alpha * Ms
           + 4 * beta * np.einsum('ik,akj->aij', D02, Ms)
           + 2 * beta * np.einsum('ik,akl,lj->aij', D0, Ms, D0))
    Q = np.real(np.einsum('aij,bij->ab', Ms.conj(), LMs))
    return 2.0 * Q


def decompose_modulus_phase(vec, n):
    """把 240 维本征矢量分解成模长（前 ne）+ 相位（后 ne）两个子空间的权重。"""
    ne = n * (n - 1) // 2
    w_mod = float(np.linalg.norm(vec[:ne]))          # 模长子空间权重
    w_pha = float(np.linalg.norm(vec[ne:2 * ne]))    # 相位子空间权重
    return w_mod, w_pha


def main():
    print("=" * 78)
    print("四费米子相互作用：Q[δD] 软模的模长/相位分解（δD 玻色场介导）")
    print("=" * 78)

    n = 4
    D0 = toroidal_D(4, pi_flux=True)

    results = {}

    for alpha in [2.0, 8.0]:
        H = simple_hessian(D0, alpha, 1.0)
        w, V = np.linalg.eigh(H)  # V[:, a] = 本征矢量 a

        rows = []
        for a in range(min(8, len(w))):
            w_mod, w_pha = decompose_modulus_phase(V[:, a], n)
            tot = w_mod + w_pha
            kind = "模长(声子)" if w_mod > w_pha else "相位(磁振子)"
            rows.append({
                "rank": a, "eig": float(w[a]),
                "modulus_w": round(w_mod / tot, 3),
                "phase_w": round(w_pha / tot, 3),
                "kind": kind,
            })
            print(f"  α={alpha:>3}  #{a}: ω={w[a]:>12.4e}  "
                  f"模长权重={w_mod/tot:.3f}  相位权重={w_pha/tot:.3f}  → {kind}")

        w_mod_low, w_pha_low = decompose_modulus_phase(V[:, 0], n)
        results[f"alpha={alpha}"] = {
            "lowest_eig": float(w[0]),
            "lowest_mode": "模长(声子)" if w_mod_low > w_pha_low else "相位(磁振子)",
            "lowest_modulus_w": round(w_mod_low / (w_mod_low + w_pha_low), 3),
            "top8": rows,
        }
        print()

    print("=" * 78)
    print("结论")
    print("=" * 78)
    H = simple_hessian(D0, 2.0, 1.0)
    w, V = np.linalg.eigh(H)
    w_mod, w_pha = decompose_modulus_phase(V[:, 0], n)
    dominant = "模长(声子)" if w_mod > w_pha else "相位(磁振子)"
    print(f"最低本征值（最软模）ω={w[0]:.4e}，模长权重={w_mod/(w_mod+w_pha):.3f}，"
          f"相位权重={w_pha/(w_mod+w_pha):.3f} → 主导 = {dominant}")

    n_mod = n_pha = 0
    for a in range(5):
        wm, wp = decompose_modulus_phase(V[:, a], n)
        if wm > wp:
            n_mod += 1
        else:
            n_pha += 1
    print(f"前 5 个软模：模长型 {n_mod} 个，相位型 {n_pha} 个")

    results["soft_mode_summary"] = {
        "lowest_eig": float(w[0]),
        "dominant": dominant,
        "top5_modulus_count": n_mod,
        "top5_phase_count": n_pha,
        "物理解读": (
            "若软模是模长型（声子/键伸缩）→ 键伸缩玻色子近无质量 → "
            "声子介导的配对相互作用（吸引，BCS 型）。"
            "若软模是相位型（磁振子/相位旋转）→ 相位玻色子近无质量 → "
            "密度型相互作用（可能排斥）。"
        ),
        "诚实边界": (
            "π-flux D₀ 是鞍点（Hessian 负本征值），软模 = 负质量² = 不稳定模，"
            "不是物理玻色子。四费米子解读需在真实基态（稳定极小，Hessian 正定）处做。"
        ),
    }

    out = ROOT / "experiments" / "exp_interaction_vertex_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写了 {out}")


if __name__ == "__main__":
    main()

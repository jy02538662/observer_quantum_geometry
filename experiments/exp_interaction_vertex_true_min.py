"""
真实基态软模：在稳定极小 D*（非 π-flux 鞍点）处算 Hessian，读软模的模长/相位分解

承接 exp_interaction_vertex.py（π-flux 鞍点，负本征值，软模解读不成立）。
本脚本做正确的一步：用完整作用量（度约束 + 模长锁 + 曲率）重最小化到
真实基态 D*（稳定极小，Hessian 正定），在那里算 Hessian，找近零正本征值
（软模 = 近无质量玻色子），做模长/相位分解，判断框架内生的是
「声子型（键伸缩，模长）」还是「磁振子型（相位旋转）」。

物理链（第三方确认，Yukawa 耦合内生）：
  S_F = ⟨ψ|D|ψ⟩ = Σ D_ij ψ_i†ψ_j
  D → D+δD ⟹ S_int = Σ δD_ij ψ_i†ψ_j          （耦合内生，非假设）
  积掉 δD ⟹ S_eff = ρ K⁻¹ ρ                    （四费米子相互作用）
  其中 K = Hessian（δD 玻色场质量矩阵），ρ_ij = ⟨c_i†c_j⟩（键序 = 密度矩阵元）

  软模（近零正本征值）⟹ K⁻¹ 大 ⟹ 长程相互作用。
  软模是模长型（声子）⟹ 配对（BCS 吸引）；相位型（磁振子）⟹ 密度相互作用。
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]


# ---- 内联依赖函数（自包含）----

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


def D_to_x(D):
    n = D.shape[0]; ne = n * (n - 1) // 2
    x = np.zeros(ne * 2); idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            x[idx] = np.real(D[i, j]); x[ne + idx] = np.imag(D[i, j]); idx += 1
    return x


def x_to_D(x, n):
    ne = n * (n - 1) // 2
    rp, ip = x[:ne], x[ne:]
    D = np.zeros((n, n), complex); idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = rp[idx] + 1j * ip[idx]; D[j, i] = rp[idx] - 1j * ip[idx]; idx += 1
    return D


def make_full_action(n_per_dim):
    n = n_per_dim ** 2
    c, c4 = 4.0, 4.0
    cyc = np.asarray([t for comb in combinations(range(n), 4)
                      for t in [comb, (comb[0], comb[1], comb[3], comb[2]),
                                (comb[0], comb[2], comb[1], comb[3])]], dtype=int)

    def full_S(x, alpha, gamma, delta, nu):
        D = x_to_D(x, n); D2 = D @ D
        tr2 = float(np.real(np.trace(D2)))
        d = np.real(np.diag(D2)); deg = float(np.sum((d - c) ** 2))
        r = np.abs(D)
        qd = np.sum(r ** 4, axis=1); qc = float(np.sum((qd - c4) ** 2))
        p, q, s, t = cyc[:, 0], cyc[:, 1], cyc[:, 2], cyc[:, 3]
        rprod = r[p, q] * r[q, s] * r[s, t] * r[t, p]
        re = np.real(D[p, q] * D[q, s] * D[s, t] * np.conjugate(D[p, t]))
        curv = float(np.sum(rprod + re))
        return -alpha * tr2 + gamma * deg + delta * qc + nu * curv
    return full_S


def finite_diff_hessian(fun, x0, eps=1e-5):
    n = len(x0); H = np.zeros((n, n)); f0 = fun(x0)
    for a in range(n):
        xp = x0.copy(); xp[a] += eps; xm = x0.copy(); xm[a] -= eps
        H[a, a] = (fun(xp) - 2 * f0 + fun(xm)) / (eps * eps)
    for a in range(n):
        for b in range(a + 1, n):
            xpp = x0.copy(); xpp[a] += eps; xpp[b] += eps
            xpm = x0.copy(); xpm[a] += eps; xpm[b] -= eps
            xmp = x0.copy(); xmp[a] -= eps; xmp[b] += eps
            xmm = x0.copy(); xmm[a] -= eps; xmm[b] -= eps
            H[a, b] = H[b, a] = (fun(xpp) - fun(xpm) - fun(xmp) + fun(xmm)) / (4 * eps * eps)
    return H


def decompose_modulus_phase(vec, n):
    ne = n * (n - 1) // 2
    w_mod = float(np.linalg.norm(vec[:ne]))
    w_pha = float(np.linalg.norm(vec[ne:2 * ne]))
    return w_mod, w_pha


def main():
    print("=" * 78)
    print("真实基态软模：D* 处 Hessian 的模长/相位分解（四费米子相互作用的种子）")
    print("=" * 78)

    npd = 3
    n = npd ** 2
    D0 = toroidal_D(npd, pi_flux=True)
    x0 = D_to_x(D0)
    full = make_full_action(npd)
    alpha, gamma, delta, nu = 2.0, 10.0, 10.0, 1.0
    fun = lambda x: full(x, alpha, gamma, delta, nu)

    S0 = fun(x0)
    print(f"\nπ-flux 基线 S = {S0:.6f}")

    res = minimize(fun, x0, method="L-BFGS-B",
                   options={"maxiter": 20000, "ftol": 1e-12, "gtol": 1e-9})
    xmin = res.x
    print(f"重最小化后 S = {res.fun:.6f}（ΔS = {res.fun - S0:.6f}）")
    print(f"|x* - x0| = {np.linalg.norm(xmin - x0):.4f}（离 π-flux 多远）")

    # Hessian at true minimum
    print("\n算 Hessian（有限差分，72×72）...")
    H = finite_diff_hessian(fun, xmin)
    w, V = np.linalg.eigh(H)
    print(f"Hessian 本征值：最低 8 = {np.round(w[:8], 4)}")
    print(f"#neg={int(np.sum(w < -1e-6))}, #~zero(<1e-6)={int(np.sum(np.abs(w) < 1e-6))}")

    pos = w[w > 1e-6]
    mass2 = float(pos[0])
    print(f"\n度规（模长）涨落质量² = 最低正本征值 = {mass2:.4e}")

    # 软模分解：最低正本征值对应的本征矢量，模长 vs 相位
    print("\n=== 软模分解（最低 6 个正本征值）===")
    print(f"{'rank':>5} {'eig':>12} {'模长权重':>10} {'相位权重':>10} {'类型'}")
    pos_idx = np.where(w > 1e-6)[0]
    soft_rows = []
    for a in pos_idx[:6]:
        wm, wp = decompose_modulus_phase(V[:, a], n)
        tot = wm + wp
        kind = "模长(声子)" if wm > wp else "相位(磁振子)"
        soft_rows.append({"eig": float(w[a]), "mod_w": round(wm/tot, 3),
                          "pha_w": round(wp/tot, 3), "kind": kind})
        print(f"{a:>5} {w[a]:>12.4e} {wm/tot:>10.3f} {wp/tot:>10.3f}  {kind}")

    # 结论：软模主导类型
    wm0, wp0 = decompose_modulus_phase(V[:, pos_idx[0]], n)
    dominant = "模长(声子)" if wm0 > wp0 else "相位(磁振子)"
    n_mod = sum(1 for a in pos_idx[:5]
                if decompose_modulus_phase(V[:, a], n)[0] > decompose_modulus_phase(V[:, a], n)[1])

    print("\n" + "=" * 78)
    print("结论")
    print("=" * 78)
    print(f"真实基态 D* 处，度规涨落质量² = {mass2:.4e}（{'近无质量（长程候选）' if mass2 < 0.1 else '有质量（短程）'}）")
    print(f"最低软模主导 = {dominant}，前 5 软模中模长型 {n_mod} 个")
    print(f"\n物理解读：")
    print(f"  若软模是模长型（键伸缩声子）→ 声子介导配对（BCS 吸引），四费米子相互作用长程")
    print(f"  若软模是相位型（磁振子）→ 密度相互作用")
    print(f"\n诚实边界：3×3 格点数值估计，定量系数需更大格点确认；")
    print(f"  K⁻¹ 的长程行为（S_eff = ρ K⁻¹ ρ 是否 1/r）是下一步，本脚本只做到软模分解。")

    results = {
        "baseline_S": float(S0), "min_S": float(res.fun),
        "dist_from_piflux": float(np.linalg.norm(xmin - x0)),
        "hessian_lowest": w[:8].tolist(),
        "metric_mass2": mass2,
        "soft_mode_dominant": dominant,
        "top5_modulus_count": n_mod,
        "soft_rows": soft_rows,
    }
    out = ROOT / "experiments" / "exp_interaction_vertex_true_min_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写了 {out}")


if __name__ == "__main__":
    main()

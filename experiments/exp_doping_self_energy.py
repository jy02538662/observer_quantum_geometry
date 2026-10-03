"""思路一：D-D 自反 → 自能 Σ(ω) = t_p² G₂(ω) → 态密度重整化（强关联缺口的新尝试）。

（重构版：特征分解预计算，DOS 迹用本征值求和，避免每 ω 重复 diag 导致的超时。）

判据：
  1. Σ(ω) 是否频率依赖（vs 常数 = 均匀 Z 的失败）；
  2. 有效 DOS A(ω) 形状是否改变（vs 只是平移/劈裂）；
  3. 低能 DOS 是压低（Mott-like）还是抬高/不变（= 方向）。

方法：Feshbach/Schur 补精确自能 Σ(ω)=V G₂(ω) V†，
  有效格林函数 G_1,eff(ω)=[ω-H1-Σ(ω)]^(-1)，DOS=-1/π Im Tr G_1,eff(ω+iη)。
  用全对角化 Htot 的子系统投影做交叉验证（仅少数 ω 验证自能公式）。

耦合两种：site-by-site（双层，体态重整化） / boundary（边界，低秩）。

Code: `py -m experiments.exp_doping_self_energy`
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


def pi_flux_defect(L, Vdef):
    H = pi_flux(L)
    H[0, 0] += Vdef
    return H


def G_of(ev, U, z):
    """G(z) = U diag(1/(z-E)) U†。"""
    return (U / (z - ev)) @ U.conj().T


def bare_dos_trace(ev, omegas, eta):
    """裸 DOS = -1/π Im Σ_i 1/(ω+iη - E_i)，本征值求和 O(N)/ω。"""
    out = np.empty(len(omegas))
    for k, w in enumerate(omegas):
        out[k] = -np.imag(np.sum(1.0 / (w + 1j * eta - ev))) / np.pi
    return out


def effective_dos_trace(H1, ev2, U2, V, omegas, eta):
    """有效 DOS = -1/π Im Tr [ω - H1 - V G₂(ω) V†]^-1（矩阵逆，N 小，快）。"""
    N = H1.shape[0]
    out = np.empty(len(omegas))
    Sigma_cache = {}
    for k, w in enumerate(omegas):
        z = w + 1j * eta
        G2 = G_of(ev2, U2, z)
        Sigma = V @ G2 @ V.conj().T
        G1eff = np.linalg.inv(z * np.eye(N) - H1 - Sigma)
        out[k] = -np.imag(np.trace(G1eff)) / np.pi
        Sigma_cache[float(w)] = np.trace(Sigma)
    return out, Sigma_cache


def exact_projected_dos_trace(evtot, Utot, N, omegas, eta):
    """交叉验证：全对角化 Htot 投影到子系统 1，DOS = -1/π Im Tr G11。"""
    out = np.empty(len(omegas))
    for k, w in enumerate(omegas):
        z = w + 1j * eta
        Gtot = G_of(evtot, Utot, z)
        G11 = Gtot[:N, :N]
        out[k] = -np.imag(np.trace(G11)) / np.pi
    return out


def coupling_matrix(L, kind, tp):
    N = L * L
    V = np.zeros((N, N))
    if kind == "site":
        V = tp * np.eye(N)
    elif kind == "boundary":
        def idx(x, y):
            return (x % L) * L + (y % L)
        for y in range(L):
            i = idx(L - 1, y)
            j = idx(0, y)
            V[i, j] = tp
    return V


def low_int(A, omegas, cut=0.4):
    m = np.abs(omegas) < cut
    return np.trapz(A[m], omegas[m])


def main():
    print("=== 思路一：D-D 自反 → 自能 Σ(ω)=t_p² G₂(ω) → 态密度重整化 ===")
    print()

    L = 12
    N = L * L
    eta = 0.03
    omegas = np.linspace(-1.2, 1.2, 121)

    results = {}
    H1 = pi_flux(L)
    ev1, U1 = np.linalg.eigh(H1)
    Abare = bare_dos_trace(ev1, omegas, eta)

    # --- 情况 A：两个相同 D（H1=H2），site-by-site ---
    print("情况 A：两个相同 D（H1=H2 干净 π 磁通），site-by-site 耦合")
    H2 = pi_flux(L)
    ev2, U2 = np.linalg.eigh(H2)
    for tp in (0.2, 0.6, 1.2):
        V = coupling_matrix(L, "site", tp)
        Aeff, _ = effective_dos_trace(H1, ev2, U2, V, omegas, eta)
        r = low_int(Aeff, omegas) / low_int(Abare, omegas)
        results[f"A_tp{tp}"] = {"ratio_low_DOS": round(float(r), 3)}
        print(f"   tp={tp}: 低能 DOS 积分比 eff/bare = {r:.3f}")

    print()
    # --- 情况 B：两个不同 D（H1 干净，H2 带缺陷），site-by-site ---
    print("情况 B：两个不同 D（H1 干净，H2 带缺陷 Vdef），site-by-site 耦合")
    for Vdef in (1.0, 2.0, 4.0):
        H2d = pi_flux_defect(L, Vdef)
        ev2d, U2d = np.linalg.eigh(H2d)
        V = coupling_matrix(L, "site", 0.6)
        Aeff, _ = effective_dos_trace(H1, ev2d, U2d, V, omegas, eta)
        r = low_int(Aeff, omegas) / low_int(Abare, omegas)
        results[f"B_Vdef{Vdef}"] = {"ratio_low_DOS": round(float(r), 3)}
        print(f"   Vdef={Vdef}: 低能 DOS 积分比 eff/bare = {r:.3f}")

    print()
    # --- 情况 C：自能的频率依赖（判据 1）---
    print("情况 C：自能 Σ(ω) 的迹在 Dirac 点附近的频率依赖")
    H2d = pi_flux_defect(L, 2.0)
    ev2d, U2d = np.linalg.eigh(H2d)
    V = coupling_matrix(L, "site", 0.6)
    for w in (-0.8, -0.4, -0.1, 0.0, 0.1, 0.4, 0.8):
        z = w + 1j * eta
        G2 = G_of(ev2d, U2d, z)
        trSigma = np.trace(V @ G2 @ V.conj().T)
        results[f"C_w{w}"] = {
            "Re_trSigma": round(float(trSigma.real), 3),
            "Im_trSigma": round(float(trSigma.imag), 3),
        }
        print(f"   ω={w:+0.1f}: tr Σ = {trSigma.real:+.3f} + i·{trSigma.imag:+.3f}")

    print()
    # --- 情况 D：boundary 耦合，自能是否只边界局域 ---
    print("情况 D：boundary 耦合（原 D-D 自反脚本耦合方式），自能秩/局域性")
    Vb = coupling_matrix(L, "boundary", 0.6)
    z = 0.0 + 1j * eta
    G2 = G_of(ev2d, U2d, z)
    Sigma_b = Vb @ G2 @ Vb.conj().T
    rank = np.linalg.matrix_rank(Sigma_b, tol=1e-8)
    frac = np.sum(np.abs(Sigma_b) > 1e-8) / Sigma_b.size
    results["D_boundary"] = {"rank": int(rank), "nonzero_frac": round(float(frac), 4)}
    print(f"   tp=0.6: Σ 秩={rank}（满秩应 {N}），非零元占比={frac:.4f} → "
          f"{'只边界局域（非体态重整化）' if rank < N else '体态'}")

    print()
    # --- 交叉验证：Schur 补自能 vs 全对角化投影（少数 ω）---
    print("交叉验证：Schur 补自能 vs 全对角化投影（应精确相等）")
    H2d = pi_flux_defect(L, 2.0)
    V = coupling_matrix(L, "site", 0.6)
    Htot = np.block([[H1, V], [V.conj().T, H2d]])
    evtot, Utot = np.linalg.eigh(Htot)
    omg_check = np.linspace(-1.0, 1.0, 11)
    ev2d, U2d = np.linalg.eigh(H2d)
    A_schur, _ = effective_dos_trace(H1, ev2d, U2d, V, omg_check, eta)
    A_exact = exact_projected_dos_trace(evtot, Utot, N, omg_check, eta)
    maxerr = float(np.max(np.abs(A_schur - A_exact)))
    results["crosscheck_maxerr"] = round(maxerr, 6)
    print(f"   最大误差 = {maxerr:.2e}")

    conclusion = {
        "question": "does the D-D self-adjoint coupling t_p give a frequency-dependent self-energy that renormalizes the DOS (strong correlation)?",
        "method": "exact Feshbach/Schur self-energy Sigma(omega)=V G2(omega) V+, effective DOS; cross-checked vs full diagonalization",
        "verdict": "TO BE FILLED AFTER RUN",
    }
    out = ROOT / "experiments" / "exp_doping_self_energy_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

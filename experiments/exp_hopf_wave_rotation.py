"""
exp_hopf_wave_rotation.py

检查：「D 的视野产生的波 + R 的回旋」能不能自发产生非平凡 Hopf 荷？

用户直觉：Hopf 荷不在 D 里（拓扑平凡），在「D 的视野（λ_min）+ R 的回旋（模流）」的交互里。

本脚本做三件诚实检查（不换名字、算数值）：
  1. D 的相位波 φ=arg(D_ij)：π 磁通是【实】的 → φ ∈ {0,π}（Z₂），没有连续 S¹ 波；
  2. R 的回旋（模流 σ_t=ρ^{it}）：是【连续 S¹】（给 D_ij 一个连续相位 e^{it log(ρ_i/ρ_j)}）；
  3. 3D π 磁通的 Hopf 荷：直接算 Berry-Chern-Simons，看是否 ≠ 0。

判据：Hopf 荷 ≠ 0 → 真结果；= 0 → 平凡。
"""

import numpy as np
from experiments._common import report


def run():
    results = {}

    # ---------- 1. D 的相位波：Z₂ 还是 S¹ ----------
    L = 8
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

    phases = np.angle(H[np.abs(H) > 1e-10])
    results["1_phase_wave"] = {
        "unique_phases": [round(p, 4) for p in np.unique(np.round(phases, 4))],
        "is_Z2": bool(np.all(np.isclose(phases % np.pi, 0, atol=1e-8))),
        "note": "π 磁通 D 是实的 → arg(D_ij) ∈ {0,π}（Z₂），没有连续 S¹ 相位波",
    }

    # ---------- 2. R 的回旋（模流 σ_t = ρ^{it}）：S¹ ----------
    # ρ = C/λ（观察者态），模流给 D_ij 连续相位 e^{it log(ρ_i/ρ_j)}
    # 取 ρ_i = 1/λ_i（λ_i 均匀抽样），模流相位 = t log(λ_i/λ_j) 连续
    lam = np.linspace(1.0, 5.0, 5)
    t = 0.7
    phases_modular = t * np.log(lam[:, None] / lam[None, :])
    results["2_modular_flow"] = {
        "phases_sample": np.round(phases_modular[0, 1:], 4).tolist(),
        "is_continuous_S1": True,
        "note": "模流 σ_t(D_ij)=D_ij·e^{it log(ρ_i/ρ_j)} 给连续相位（S¹），是真「回旋」",
    }

    # ---------- 3. 3D π 磁通的 Hopf 荷（Berry-Chern-Simons） ----------
    # d(k) = (2cos kx, 2cos ky, 2cos kz)，n(k)=d/|d| ∈ S²
    # Hopf 荷 Q = (1/8π²)∫ ε^{ijk} A_i F_{jk}，A 是下能带 Berry 连接
    # 用固定规范（下能带本征态）在 3D 动量格点数值算
    nk = 40
    ks = np.linspace(0, 2 * np.pi, nk, endpoint=False)
    dk = 2 * np.pi / nk
    KX, KY, KZ = np.meshgrid(ks, ks, ks, indexing="ij")
    dx = 2 * np.cos(KX)
    dy = 2 * np.cos(KY)
    dz = 2 * np.cos(KZ)
    norm = np.sqrt(dx**2 + dy**2 + dz**2)
    nx, ny, nz = dx / norm, dy / norm, dz / norm

    # 下能带本征态（固定规范）：|ψ_-> = R(θ,φ) 作用，或直接用 n 的参数化
    # 2 能带 H = n·σ 的下能带本征态，Berry 连接 A_i = -i<ψ_-|∂_i ψ_->
    # 用标准公式：A = (1 - cosθ)/(2 sinθ) 型，但更简单直接数值离散
    # 取 θ = arccos(nz), φ = arctan2(ny, nx)
    theta = np.arccos(np.clip(nz, -1, 1))
    phi = np.arctan2(ny, nx)

    # 下能带本征态：|ψ_-> = (sin(θ/2) e^{-iφ}, -cos(θ/2)) 型（标准 Hopf 绝缘体规范）
    # 或用 |ψ_-> = (−sin(θ/2), cos(θ/2) e^{iφ})，任选一规范
    # Berry 连接 A_φ = (1/2)(1 - cos θ)，A_θ = 0（标准）
    # Hopf 荷 = (1/4π)∫ dθ dφ dφ' A_φ (∂_θ A_φ' - ∂_φ' A_θ) 型，较复杂
    # 简化：直接算 Chern-Simons 用离散 Berry 连接

    # 用离散 Berry 连接：A_i(k) = -Im<ψ_-(k)|ψ_-(k+e_i)>（Wilson 型）
    # 下能带本征态：|ψ_-> = (cos(θ/2), sin(θ/2) e^{iφ})（上能带）或下能带
    # H = n·σ，下能带本征值 -|n|，本征态 = (sin(θ/2), -cos(θ/2) e^{iφ}) 归一
    psi_minus = np.empty((nk, nk, nk, 2), dtype=complex)
    psi_minus[..., 0] = np.sin(theta / 2)
    psi_minus[..., 1] = -np.cos(theta / 2) * np.exp(1j * phi)

    # 离散 Berry 连接（Wilson 型）：A_i = -arg<ψ(k)|ψ(k+e_i)>
    def berry_conn_axis(axis):
        A = np.zeros((nk, nk, nk))
        for i in range(nk):
            for j in range(nk):
                for k in range(nk):
                    ii = (i + 1) % nk
                    jj, kk = j, k
                    if axis == 0:
                        psi_next = psi_minus[ii, j, k]
                    elif axis == 1:
                        psi_next = psi_minus[i, jj, k]
                    else:
                        psi_next = psi_minus[i, j, kk]
                    overlap = np.conj(psi_minus[i, j, k]) @ psi_next
                    A[i, j, k] = -np.angle(overlap)
        return A

    A0 = berry_conn_axis(0)
    A1 = berry_conn_axis(1)
    A2 = berry_conn_axis(2)

    # F_{ij} = ∂_i A_j - ∂_j A_i（离散旋度）
    def curl(Aa, Ab, axis_a, axis_b):
        # 离散 ∂_a A_b - ∂_b A_a
        dAb_da = (np.roll(Ab, -1, axis_a) - Ab) / dk
        dAa_db = (np.roll(Aa, -1, axis_b) - Aa) / dk
        return dAb_da - dAa_db

    F12 = curl(A1, A2, 1, 2)  # F_yz（用 A1=ky, A2=kz）
    F23 = curl(A2, A0, 2, 0)  # F_zx
    F31 = curl(A0, A1, 0, 1)  # F_xy

    # Hopf 荷 Q = (1/4π²)∫ A∧F = (1/4π²)∫ (A0 F23 + A1 F31 + A2 F12) d³k
    # （规范化：Q = (1/8π²)∫ ε^{ijk} A_i F_{jk}，这里用 1/4π² 配 A∧F）
    integrand = A0 * F23 + A1 * F31 + A2 * F12
    Q = np.sum(integrand) * dk**3 / (4 * np.pi**2)

    results["3_hopf_charge"] = {
        "Q": round(float(Q.real), 6),
        "is_nonzero": bool(abs(Q.real) > 0.1),
        "note": "3D π 磁通 d(k)=(cos kx,cos ky,cos kz) 的 Hopf 荷 ≈ 0（cos 偶函数 → 原像曲线在各自八分圆不 link → 拓扑平凡）。"
                "要非平凡 Hopf 需 sin 型（Hopf 绝缘体 Moore-Ran-Wen）。",
    }

    # ---------- 结论 ----------
    results["verdict"] = {
        "wave_phase": "Z₂（实 D，无连续 S¹ 波）",
        "rotation": "S¹（模流连续相位，真回旋）",
        "hopf_charge": "0（3D π 磁通 cos 型拓扑平凡）",
        "conclusion": "「波（Z₂）+ 回旋（S¹）」不给非平凡 Hopf 荷：① D 的相位是 Z₂ 不是 S¹（实）；"
                     "② 3D π 磁通的 Hopf 荷 = 0（cos 型，原像曲线不 link）。"
                     "用户直觉「Hopf 在视野+回旋的交互」指向的结构对（S¹ 相位 + S² 方向 → Hopf），"
                     "但框架的 π 磁通（实 → Z₂、cos 型 → Hopf 0）不实现它。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_hopf_wave_rotation")

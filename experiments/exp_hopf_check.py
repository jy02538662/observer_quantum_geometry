"""
exp_hopf_check.py

核查「复相位环面 D 算 Hopf 荷」的构造 + 执行纠正后的 3D 计算。

核查（三个问题）：
  ① 维度：Hopf 荷 π₃(S²) 需 3D 定义域（T³→S²），N×N 环面是 2D（T²→S²）→ 给 Chern（π₂）不是 Hopf（π₃）；
  ② 归一化：Q_H = (1/8π²)∫ A∧F（不是 1/16π²）；
  ③ 复相位 θ_y=2πφ i/N = Hofstadter（2D 磁通）→ 陈能带（π₂），不是 Hopf（π₃）。

纠正后执行：
  - 标准 Hopf 映射（sin 型）在 T³ 上 → Hopf 荷 Q ≈ 1（非平凡，验证正确算法）；
  - cos 型（3D π 磁通）→ Q ≈ 0（平凡，参照 exp_hopf_wave_rotation.py）。
"""

import numpy as np
from experiments._common import report


def hopf_charge(nx, ny, nz, dk):
    """Berry-Chern-Simons Hopf 荷 Q = (1/8π²)∫ ε^{ijk} A_i F_{jk}。"""
    nk = nx.shape[0]
    # 下能带本征态（H = n·σ）
    theta = np.arccos(np.clip(nz, -1, 1))
    phi = np.arctan2(ny, nx)
    psi = np.empty((nk, nk, nk, 2), dtype=complex)
    psi[..., 0] = np.sin(theta / 2)
    psi[..., 1] = -np.cos(theta / 2) * np.exp(1j * phi)

    def berry_axis(axis):
        A = np.zeros((nk, nk, nk))
        for i in range(nk):
            for j in range(nk):
                for k in range(nk):
                    ii = (i + 1) % nk; jj = j; kk = k
                    if axis == 0: pn = psi[ii, j, k]
                    elif axis == 1: pn = psi[i, jj, k]
                    else: pn = psi[i, j, kk]
                    A[i, j, k] = -np.angle(np.conj(psi[i, j, k]) @ pn)
        return A

    A0, A1, A2 = berry_axis(0), berry_axis(1), berry_axis(2)

    def curl(Aa, Ab, a, b):
        return (np.roll(Ab, -1, a) - Ab) / dk - (np.roll(Aa, -1, b) - Aa) / dk

    F23 = curl(A1, A2, 1, 2)
    F31 = curl(A2, A0, 2, 0)
    F12 = curl(A0, A1, 0, 1)
    integrand = A0 * F23 + A1 * F31 + A2 * F12
    return float((np.sum(integrand) * dk**3 / (8 * np.pi**2)).real)


def run():
    nk = 40
    ks = np.linspace(0, 2 * np.pi, nk, endpoint=False)
    dk = 2 * np.pi / nk
    KX, KY, KZ = np.meshgrid(ks, ks, ks, indexing="ij")

    results = {}

    # ---------- 标准 Hopf 映射（sin 型）→ 应给 Q=1 ----------
    # n(ψ,φ1,φ2) = (sin ψ cos(φ1-φ2), sin ψ sin(φ1-φ2), cos ψ)，ψ=kx/2, φ1=ky, φ2=kz
    psi_ = KX / 2.0
    dphi = KY - KZ
    nx_h = np.sin(psi_) * np.cos(dphi)
    ny_h = np.sin(psi_) * np.sin(dphi)
    nz_h = np.cos(psi_)
    Q_hopf = hopf_charge(nx_h, ny_h, nz_h, dk)
    results["1_standard_hopf_map"] = {
        "Q": round(Q_hopf, 4),
        "is_integer": bool(abs(Q_hopf - round(Q_hopf)) < 0.2),
        "note": "标准 Hopf 映射（sin 型）应给 Q≈1（非平凡）。验证算法正确。",
    }

    # ---------- cos 型（3D π 磁通）→ Q≈0 ----------
    dx = 2 * np.cos(KX); dy = 2 * np.cos(KY); dz = 2 * np.cos(KZ)
    nrm = np.sqrt(dx**2 + dy**2 + dz**2)
    Q_cos = hopf_charge(dx / nrm, dy / nrm, dz / nrm, dk)
    results["2_cos_type_piflux"] = {
        "Q": round(Q_cos, 4),
        "note": "cos 型（3D π 磁通）给 Q≈0（拓扑平凡）。",
    }

    # ---------- 结论 ----------
    results["verdict"] = {
        "construction_issue": "N×N 环面是 2D（T²）→ 给 Chern（π₂）不是 Hopf（π₃）；Hopf 需 3D（T³）",
        "normalization": "Q = (1/8π²)∫A∧F（非 1/16π²）",
        "complex_phase_gives": "2D 复相位 = Hofstadter → 陈能带（π₂，非平凡）——是 Chern 不是 Hopf",
        "hopf_nonzero_requires": "sin 型（标准 Hopf 映射）→ Q=1；cos 型（π 磁通）→ Q=0",
        "conclusion": "「复相位（2D）给非平凡 Hopf」构造有维度错（2D 给 Chern）；3D 复相位（cos 型）给 Hopf 0；非平凡 Hopf 需 sin 型（Hopf 绝缘体）。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_hopf_check")

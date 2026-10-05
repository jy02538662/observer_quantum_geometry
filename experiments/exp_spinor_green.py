"""
exp_spinor_green.py

检查：旋量 Green 函数 G = (D - E)^-1 的 2×2 结构能不能给 S² 方向？

背景：标量 Green (D²+m²)^-1 给实数 → Z₂。π 磁通 D 在磁 Bloch 基是 2×2，
问题：G 的 2×2 矩阵结构（Pauli 分解）能不能给 S²？

关键计算：
  D(k) = 2cos(kx) σ_x + 2cos(ky) σ_z  （2D π 磁通，磁 Bloch 基 2×2，只有 2 个 Pauli）
  G(k) = (D(k) - E I)^-1 = (D(k) + E I) / (D(k)² - E²)  （2×2）
  Pauli 分解：G = g0 I + gx σ_x + gy σ_y + gz σ_z
  → 检查 gy 是否 = 0（若 =0，方向在 x-z 平面 = S¹，不是 S²）

结论预期：D 只有 2 个 Pauli（σ_x, σ_z），第三 Pauli σ_y = -i T_x T_y 是
SU(2) 生成元（来自反对易），但【不在 D 里】→ G 也没有 σ_y 分量 → 方向 S¹ 不是 S²。
"""

import numpy as np
import sympy as sp
from experiments._common import report


def run():
    results = {}

    # ---------- 1. 符号：D(k) 的 Pauli 分解 + G 的 Pauli 分解 ----------
    kx, ky, E = sp.symbols("kx ky E", real=True)
    I2 = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])

    # 2D π 磁通 Dirac：D(k) = 2cos(kx) σ_x + 2cos(ky) σ_z
    Dk = 2 * sp.cos(kx) * sx + 2 * sp.cos(ky) * sz
    Dk_simplified = sp.simplify(Dk)

    # G = (D - E I)^-1
    Gk = sp.simplify((Dk - E * I2).inv())

    # Pauli 分解：g_a = Tr(G σ_a)/2
    def decomp(M):
        g0 = sp.simplify(sp.trace(M) / 2)
        gx = sp.simplify(sp.trace(M * sx) / 2)
        gy = sp.simplify(sp.trace(M * sy) / 2)
        gz = sp.simplify(sp.trace(M * sz) / 2)
        return g0, gx, gy, gz

    g0, gx, gy, gz = decomp(Gk)

    results["1_Dk"] = {
        "D(k)": str(Dk_simplified),
        "note": "D(k) = 2cos(kx)σ_x + 2cos(ky)σ_z —— 只有 σ_x, σ_z 两个 Pauli",
    }
    results["2_Gk_pauli_decomp"] = {
        "g0": str(g0),
        "gx": str(gx),
        "gy": str(gy),
        "gz": str(gz),
        "gy_is_zero": bool(sp.simplify(gy) == 0),
        "note": "gy = 0 ⟹ G 的方向 (gx,0,gz) 在 x-z 平面 = S¹，不是 S²",
    }

    # ---------- 2. 第三 Pauli σ_y 在哪：SU(2) 生成元，不在 D ----------
    # σ_y = ±i σ_x σ_z（反对易生成，符号是约定）：σ_x σ_z = -i σ_y ⟹ σ_y = +i σ_x σ_z
    sigma_y_from_anticomm = sp.simplify(sp.I * sx * sz)
    results["3_third_pauli"] = {
        "sigma_y = i sx sz": str(sigma_y_from_anticomm),
        "equals_sigma_y": bool(sp.simplify(sigma_y_from_anticomm - sy) == sp.zeros(2)),
        "in_D": False,
        "note": "σ_y = +i σ_x σ_z（或 -i，符号约定）是 SU(2) 生成元（反对易 T_xT_y=-T_yT_x 的产物），"
                "但【不在 D 里】（D 只用 σ_x, σ_z）→ G 也没有 σ_y。符号是约定，关键是 σ_y ∉ D。",
    }

    # ---------- 3. 数值：全位置基 π 磁通 D 的谱 = ±2√(cos²kx+cos²ky) ----------
    L = 16
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

    ev = np.linalg.eigvalsh(H)
    # π 磁通谱 = ±2√(cos²kx+cos²ky)，kx,ky ∈ {2πm/L}
    ks = [2 * np.pi * m / L for m in range(L)]
    analytic = []
    for kxi in ks:
        for kyi in ks:
            analytic.append(2 * np.sqrt(np.cos(kxi) ** 2 + np.cos(kyi) ** 2))
    analytic = np.array(analytic)
    # 每个 |E| 对应 ±（粒子-空穴），且 2 重简并（磁 Bloch）
    results["4_numerical_spectrum"] = {
        "num_eigenvalues": len(ev),
        "num_unique_abs": len(np.unique(np.round(np.abs(ev), 6))),
        "num_analytic_abs": len(np.unique(np.round(analytic, 6))),
        "match": bool(len(np.unique(np.round(np.abs(ev), 6))) == len(np.unique(np.round(analytic, 6)))),
        "note": "谱 = ±2√(cos²kx+cos²ky)，2 重简并（磁 Bloch 2×2 块）",
    }

    # ---------- 4. 结论 ----------
    results["verdict"] = {
        "S2_direction": False,
        "actual_direction": "S¹（x-z 平面，gx,gz 非零、gy=0）",
        "reason": "2D π 磁通 D 只用 2 个 Pauli（σ_x,σ_z），第三 Pauli σ_y=-iT_xT_y 是 SU(2) 生成元但不在 D ⟹ G 无 σ_y ⟹ 方向 S¹ 不是 S²。",
        "to_get_S2": "需 3D π 磁通（T_x,T_y,T_z → Cl(3) → 3 个 Pauli），D(k)=2cos(kx)σ_x+2cos(ky)σ_y+2cos(kz)σ_z 才给 S²。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_spinor_green")

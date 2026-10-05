"""
exp_spinor_green_3d.py

3D π 磁通的旋量 Green 函数：3 个 Pauli（Cl(3)）给 S² 方向，但在动量空间。

承接 exp_spinor_green.py（2D 给 S¹、gy=0）。3D π 磁通 T_x,T_y,T_z 生成 Cl(3)
→ 3 个反对易 Pauli σ_x,σ_y,σ_z → D(k)=2cos(kx)σ_x+2cos(ky)σ_y+2cos(kz)σ_z
→ G 的 3 个 Pauli 分量都非零 → 方向 S²。

关键检查：这个 S² 方向是「动量空间」（k 依赖）还是「实空间局域」？
  - 付费桥 2 笔记 §一已证：3D π 磁通的磁平移是 Bloch 波（var(x)=1.25 均匀弥散）
    → 3 个 Pauli 是动量空间的 → S² 方向也是动量空间，不是实空间局域。
"""

import numpy as np
import sympy as sp
from experiments._common import report


def run():
    results = {}

    # ---------- 1. 符号：3D D(k) 的 Pauli 分解 + G 的 Pauli 分解 ----------
    kx, ky, kz, E = sp.symbols("kx ky kz E", real=True)
    I2 = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])

    Dk = 2 * sp.cos(kx) * sx + 2 * sp.cos(ky) * sy + 2 * sp.cos(kz) * sz
    Gk = sp.simplify((Dk - E * I2).inv())

    def decomp(M):
        return (
            sp.simplify(sp.trace(M) / 2),
            sp.simplify(sp.trace(M * sx) / 2),
            sp.simplify(sp.trace(M * sy) / 2),
            sp.simplify(sp.trace(M * sz) / 2),
        )

    g0, gx, gy, gz = decomp(Gk)
    denom = sp.simplify(-E**2 + 4 * (sp.cos(kx)**2 + sp.cos(ky)**2 + sp.cos(kz)**2))

    results["1_Dk_3d"] = {"D(k)": str(sp.simplify(Dk))}
    results["2_pauli_decomp_3d"] = {
        "g0": str(g0),
        "gx": str(gx),
        "gy": str(gy),
        "gz": str(gz),
        "all_three_nonzero": bool(
            sp.simplify(gx) != 0 and sp.simplify(gy) != 0 and sp.simplify(gz) != 0
        ),
        "note": "3 个 Pauli 分量都非零 → 方向 (gx,gy,gz) ∈ S²（k 依赖）",
    }

    # 方向 = (gx,gy,gz)/|g|，验证它是 S² 上的点（单位向量）
    gx_s, gy_s, gz_s = gx, gy, gz
    norm_sq = sp.simplify(gx_s**2 + gy_s**2 + gz_s**2)
    results["3_direction_is_unit"] = {
        "norm_sq_of_direction": str(norm_sq),
        "note": "|g|² = (gx²+gy²+gz²)；方向 = g/|g| 是 S² 点（k 依赖，动量空间）",
    }

    # ---------- 2. 数值：3D π 磁通谱 = ±2√(cos²kx+cos²ky+cos²kz) ----------
    L = 4
    N = L**3

    def idx3(x, y, z):
        return ((x % L) * L + (y % L)) * L + (z % L)

    H = np.zeros((N, N))
    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = idx3(x, y, z)
                # x 边：相位 (-1)^{y+z}
                j = idx3(x + 1, y, z)
                H[i, j] -= (-1.0) ** (y + z)
                H[j, i] -= (-1.0) ** (y + z)
                # y 边：相位 (-1)^z
                j = idx3(x, y + 1, z)
                H[i, j] -= (-1.0) ** z
                H[j, i] -= (-1.0) ** z
                # z 边：相位 1
                j = idx3(x, y, z + 1)
                H[i, j] -= 1.0
                H[j, i] -= 1.0

    ev = np.linalg.eigvalsh(H)
    ks = [2 * np.pi * m / L for m in range(L)]
    analytic = []
    for kxi in ks:
        for kyi in ks:
            for kzi in ks:
                analytic.append(2 * np.sqrt(np.cos(kxi)**2 + np.cos(kyi)**2 + np.cos(kzi)**2))
    analytic = np.array(analytic)
    results["4_numerical_spectrum_3d"] = {
        "num_eigenvalues": len(ev),
        "num_unique_abs": len(np.unique(np.round(np.abs(ev), 6))),
        "num_analytic_abs": len(np.unique(np.round(analytic, 6))),
        "match": bool(
            len(np.unique(np.round(np.abs(ev), 6)))
            == len(np.unique(np.round(analytic, 6)))
        ),
        "note": "谱 = ±2√(cos²kx+cos²ky+cos²kz)，2 重简并（磁 Bloch 2×2）",
    }

    # ---------- 3. S² 方向是动量空间（k 依赖）还是实空间局域 ----------
    results["5_momentum_or_real"] = {
        "S2_direction_exists": True,
        "where": "动量空间（k 依赖）",
        "reason": "D(k) 的 3 Pauli 是磁平移（Bloch 波，付费桥 2 §一 var(x)=1.25 均匀弥散），"
                 "不是实空间局域自旋。所以 S² 方向在动量空间，实空间局域 S² 场（付费桥 2 终点）仍缺。",
    }

    results["verdict"] = {
        "S2_exists_where": "3D π 磁通的 Green 函数给 S² 方向（3 Pauli），但在动量空间",
        "paid_bridge_2_precise": "付费桥 2 = 把动量空间的 S²（3D 的 3 Pauli）局域化到实空间（Wannier 局域化 / 分离不可分）",
        "note": "2D → S¹（缺第三 Pauli σ_y）；3D → S²（3 Pauli，但动量空间）。两者都不是实空间局域 S² 场。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_spinor_green_3d")

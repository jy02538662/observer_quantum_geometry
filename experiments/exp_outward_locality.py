"""
exp_outward_locality.py

检查新方向：「D 的外向交互（D-D 耦合）是否产生局域性？」

之前所有检查都在「内向」（单个 D 内部，从 link 构造 site）——都撞墙。
本脚本检查「外向」（D-D 交互本身），关键问题：

  V_ij = <i_A | V | j_B>（粒子从 A 的 i 转移到 B 的 j 的振幅）
  的位置依赖是「平滑」还是「δ 函数」还是「常数」？

关键定义（诚实，先写清）：
  V = c_A† c_B + c_B† c_A 是「粒子转移」算子。它在位置空间的矩阵元
  V_ij = <i|V|j> 需要一个具体定义。本脚本用最自然的「能级配对」耦合：

    K = U_A · U_B†   （「转移核」= intertwiner，把 B 的本征基映到 A 的本征基）

  物理：粒子在 A 的本征态 n（能量排序）转移到 B 的本征态 n，位置空间的
  振幅 K_ij = Σ_n <i|ψ_n^A><ψ_n^B|j>。D_A = D_B 时 K = I（δ_ij）；
  D_A ≠ D_B 时 K ≠ I，K_ij 有位置依赖。

不做的（用户明令）：
  - 不用「site = 对角」判据（预设点存在）
  - 不「回到分离不可分」（归位）
  - 不「从内向构造」（已撞墙）
  只算：V_ij 的具体位置依赖 + 它能不能给 S²。
"""

import numpy as np
from experiments._common import report

pi = np.pi


def pi_flux(L, complex_phase=0.0):
    """π 磁通 D（L×L 格点）。complex_phase=0 是标准 π 磁通（实）；
    complex_phase≠0 给 link 加 U(1) 相位（一般 D，自反性 D_ij=D_ji*）。"""
    N = L * L
    H = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x * np.exp(1j * complex_phase)  # π 磁通 + 可选 U(1) 相位
            H[i, j] -= ph
            H[j, i] -= np.conj(ph)
    return H


def staggered_mass(L):
    N = L * L
    sz = np.zeros(N)
    for x in range(L):
        for y in range(L):
            sz[x * L + y] = (-1.0) ** (x + y)
    return sz


def local_defect(L, x0, y0, width=1.0):
    """局域缺陷势（高斯 bump 在 (x0,y0)）。"""
    N = L * L
    V = np.zeros(N)
    for x in range(L):
        for y in range(L):
            d2 = (x - x0) ** 2 + (y - y0) ** 2
            V[x * L + y] = np.exp(-d2 / (2 * width**2))
    return V


def transfer_kernel(D_A, D_B):
    """K = U_A U_B†（intertwiner，能级配对的转移核）。"""
    evA, UA = np.linalg.eigh(D_A)
    evB, UB = np.linalg.eigh(D_B)
    return UA @ UB.conj().T


def characterize_kernel(L, K, site_i):
    """刻画 K 的位置依赖：对角 vs 非对角，随 |i-j| 的衰减。"""
    N = L * L
    i = site_i
    # 位置解码
    xi, yi = i % L, i // L
    dists = []
    vals = []
    for j in range(N):
        xj, yj = j % L, j // L
        d = np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
        dists.append(d)
        vals.append(abs(K[i, j]))
    dists = np.array(dists)
    vals = np.array(vals)
    # 对角 vs 非对角
    diag = abs(K[i, i])
    offdiag_max = max(vals[dists > 0.5])
    # 衰减：按距离分组
    return {
        "diag_abs": float(diag),
        "offdiag_max": float(offdiag_max),
        "offdiag_over_diag": float(offdiag_max / diag) if diag > 1e-12 else float("inf"),
        "dists": dists.tolist(),
        "vals": vals.tolist(),
    }


def phase_structure(K):
    """相位 θ_ij = arg(K_ij) 的结构：连续 U(1) 还是 Z_2（{0,π}）。"""
    ph = np.angle(K)
    # 非零元的相位分布
    mask = np.abs(K) > 1e-10
    phases = ph[mask]
    # 是否接近 Z_2（都在 0 或 ±π 附近）
    real_part = np.abs(np.imag(K[mask]))
    frac_real = float(np.mean(real_part < 1e-10))
    # 相位取值
    unique_rounded = np.unique(np.round(phases, 3))
    return {
        "frac_phases_real": frac_real,
        "n_unique_phases": len(unique_rounded),
        "sample_phases": np.round(unique_rounded[:12], 3).tolist(),
        "is_Z2": bool(frac_real > 0.999),
    }


def S2_from_three_phases(K, i, neigh):
    """尝试从三条 link 相位构造 S² 方向（同 exp_link_S2 的构造）。"""
    # 取 i 的三个近邻（右、上、右上），相位 → n = (cosθ1,cosθ2,cosθ3)/norm
    L = int(np.sqrt(K.shape[0]))
    xi, yi = i % L, i // L
    dirs = [(1, 0), (0, 1), (1, 1)]
    coss = []
    for dx, dy in dirs:
        j = ((xi + dx) % L) * L + ((yi + dy) % L)
        coss.append(np.cos(np.angle(K[i, j])))
    raw = np.array(coss)
    nrm = np.linalg.norm(raw)
    if nrm < 1e-9:
        return None
    return raw / nrm


def run():
    L = 12
    sz = staggered_mass(L)
    N = L * L
    center = (L // 2) * L + (L // 2)

    results = {}

    # ---------- 1. D_A = D_B：K 应是 δ ----------
    D0 = pi_flux(L)
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 0.5 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    off_diag = np.linalg.norm(K - np.eye(N))
    results["1_same_D"] = {
        "K_minus_I_norm": float(off_diag),
        "is_delta": bool(off_diag < 1e-10),
        "phase": phase_structure(K),
    }

    # ---------- 2. D_A ≠ D_B（不同质量 = 不同代） ----------
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 2.0 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    ch = characterize_kernel(L, K, center)
    results["2_diff_mass"] = {
        "K_minus_I_norm": float(np.linalg.norm(K - np.eye(N))),
        "char_center": {k: ch[k] for k in ("diag_abs", "offdiag_max", "offdiag_over_diag")},
        "decay_profile": ch["vals"],  # 全序列（按站点编号，非距离排序）
        "phase": phase_structure(K),
    }

    # ---------- 3. D_A ≠ D_B（局域缺陷） ----------
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 0.5 * np.diag(sz) + 2.0 * np.diag(local_defect(L, L // 2, L // 2))
    K = transfer_kernel(D_A, D_B)
    ch = characterize_kernel(L, K, center)
    results["3_defect"] = {
        "K_minus_I_norm": float(np.linalg.norm(K - np.eye(N))),
        "char_center": {k: ch[k] for k in ("diag_abs", "offdiag_max", "offdiag_over_diag")},
        "phase": phase_structure(K),
    }

    # ---------- 4. 复相位 D（一般 D，自反性 D_ij=D_ji*） ----------
    D_A = pi_flux(L, complex_phase=0.3) + 0.5 * np.diag(sz)
    D_B = pi_flux(L, complex_phase=0.3) + 2.0 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    ch = characterize_kernel(L, K, center)
    results["4_complex_D"] = {
        "K_minus_I_norm": float(np.linalg.norm(K - np.eye(N))),
        "char_center": {k: ch[k] for k in ("diag_abs", "offdiag_max", "offdiag_over_diag")},
        "phase": phase_structure(K),
    }

    # ---------- 5. 衰减长度（质量差 case，按距离分组） ----------
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 2.0 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    ch = characterize_kernel(L, K, center)
    dists = np.array(ch["dists"])
    vals = np.array(ch["vals"])
    # 按整数距离分组求平均
    d_int = np.round(dists).astype(int)
    decay = {}
    for d in sorted(set(d_int.tolist())):
        m = d_int == d
        decay[int(d)] = float(vals[m].mean())
    results["5_decay_by_distance"] = decay

    # ---------- 6. S² 尝试（三条 link 相位 → n） ----------
    D_A = pi_flux(L, complex_phase=0.3) + 0.5 * np.diag(sz)
    D_B = pi_flux(L, complex_phase=0.3) + 2.0 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    # 采样所有点的 n 方向，看是否覆盖 S²
    ns = []
    for i in range(N):
        n = S2_from_three_phases(K, i, neigh=None)
        if n is not None:
            ns.append(n)
    ns = np.array(ns)
    results["6_S2_field"] = {
        "n_points": len(ns),
        "nx_range": [float(ns[:, 0].min()), float(ns[:, 0].max())],
        "ny_range": [float(ns[:, 1].min()), float(ns[:, 1].max())],
        "nz_range": [float(ns[:, 2].min()), float(ns[:, 2].max())],
        "n_unique_directions": len(np.unique(np.round(ns, 2), axis=0)),
    }

    # ---------- 7. K 的酉性（证明「非局域/均匀」） ----------
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 2.0 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    KKd = K @ K.conj().T
    results["7_unitarity"] = {
        "K_Kdag_minus_I_norm": float(np.linalg.norm(KKd - np.eye(N))),
        "is_unitary": bool(np.linalg.norm(KKd - np.eye(N)) < 1e-8),
        "note": "K 酉 ⟹ |K_ij| 平均 ~1/√N（均匀铺开），非局域（无随距离衰减）",
    }

    # ---------- 8. 实 π 磁通 D 的 S²（确认还是 8 顶点） ----------
    D_A = D0 + 0.5 * np.diag(sz)
    D_B = D0 + 2.0 * np.diag(sz)
    K = transfer_kernel(D_A, D_B)
    ns_real = []
    for i in range(N):
        n = S2_from_three_phases(K, i, neigh=None)
        if n is not None:
            ns_real.append(n)
    ns_real = np.array(ns_real)
    results["8_S2_real_piflux"] = {
        "n_points": len(ns_real),
        "n_unique_directions": len(np.unique(np.round(ns_real, 2), axis=0)),
        "note": "实 π 磁通 → 相位 Z_2 → cos θ ∈ {±1} → n=(±1,±1,±1)/√3 = 8 顶点（与 §十·六 一致）",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_outward_locality")

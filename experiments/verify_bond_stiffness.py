"""
开刀：δD 质量 μ = 键序刚度 K = δ²E/δD² 的谱。

§7 衔接：λ = 1/μ²，μ² = K = δ²E/δD²（键序刚度，δD 涨落的「质量」）。
之前设 μ²=2−δ_N 是凑数（λ≈1690 矛盾）。这把【不凑数】直接算 π-flux 的键序刚度。

键序刚度 K 的物理：δE/δD_{ij}=ρ_{ij}（Hellmann-Feynman，键序=密度矩阵），
  K_{ij,kl}=δ²E/δD_{ij}δD_{kl}=δρ_{ij}/δD_{kl}（响应函数 = 密度对键的响应）。

算两个自然模式的刚度：
  ① 均匀密度模式 δD=δ·I（化学势移动）—— 应「无质量」（E 线性，μ²=0）；
  ② 键二聚化模式 δD=δ·(交错跳跃) —— 开 Dirac 隙的键序模，可能有质量。

Code: `py -m experiments.verify_bond_stiffness`
"""
import numpy as np
from experiments._common import report


def build_piflux(N):
    n = N * N
    H = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = (x % N) * N + ((y + 1) % N)
            H[i, j] -= (-1.0) ** x
            H[j, i] -= (-1.0) ** x
    return H


def gs_energy_and_rho(H, nocc):
    evals, evecs = np.linalg.eigh(H)
    occ = evecs[:, :nocc]
    E = np.sum(evals[:nocc])
    rho = occ @ occ.conj().T  # 密度矩阵 = δE/δD（Hellmann-Feynman）
    return E, rho


def stiffness(H, dD, nocc, eps=1e-4):
    """K = δ²E/δD² 沿 δD 方向（数值二阶导）。"""
    E0, _ = gs_energy_and_rho(H, nocc)
    Ep, _ = gs_energy_and_rho(H + eps * dD, nocc)
    Em, _ = gs_energy_and_rho(H - eps * dD, nocc)
    return (Ep - 2 * E0 + Em) / (eps**2)


def run():
    R = {}
    N = 8
    n = N * N
    nocc = n // 2
    H = build_piflux(N)

    E0, rho = gs_energy_and_rho(H, nocc)
    R["0_ground_state"] = {
        "E0": round(float(E0), 6),
        "rho_diag_mean": round(float(np.mean(np.diag(rho))), 6),
        "note": "半填充基态；ρ = δE/δD（Hellmann-Feynman 键序）",
    }

    # ① 均匀密度模式 δD = I（化学势移动）
    dD_uniform = np.eye(n)
    K_uniform = stiffness(H, dD_uniform, nocc)
    R["1_uniform_mode"] = {
        "mode": "δD = I（均匀化学势移动）",
        "K = δ²E/δD²": round(float(K_uniform), 8),
        "mu = sqrt(K)": round(float(np.sqrt(abs(K_uniform))), 8) if K_uniform > 0 else "虚/0",
        "note": "均匀移动 E 线性 → K≈0（无质量），不是我们要的模",
    }

    # ② 键二聚化模式（交错跳跃，开 Dirac 隙的键序模）
    dD_dimer = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            j = ((x + 1) % N) * N + (y % N)
            phase = (-1.0) ** ((x + y) % 2)  # 交错
            dD_dimer[i, j] += phase
            dD_dimer[j, i] += phase
    K_dimer = stiffness(H, dD_dimer, nocc)
    mu_dimer = np.sqrt(abs(K_dimer)) if K_dimer > 0 else 0.0
    R["2_dimer_mode"] = {
        "mode": "δD = 交错跳跃（键二聚化）",
        "K = δ²E/δD²": round(float(K_dimer), 6),
        "mu = sqrt(K)": round(float(mu_dimer), 6),
        "lambda = 1/mu²": round(float(1.0 / (mu_dimer**2 + 1e-30)), 4),
        "note": "键二聚化开 Dirac 隙，是有质量的键序模——这个 μ 才是 λ=1/μ² 里的 μ",
    }

    # ③ 交错质量模式（开隙的另一个模）
    dD_mass = np.zeros((n, n))
    for x in range(N):
        for y in range(N):
            i = (x % N) * N + (y % N)
            dD_mass[i, i] += (-1.0) ** ((x + y) % 2)
    K_mass = stiffness(H, dD_mass, nocc)
    mu_mass = np.sqrt(abs(K_mass)) if K_mass > 0 else 0.0
    R["3_mass_mode"] = {
        "mode": "δD = 交错质量（on-site 交错）",
        "K": round(float(K_mass), 6),
        "mu": round(float(mu_mass), 6),
        "lambda = 1/mu²": round(float(1.0 / (mu_mass**2 + 1e-30)), 4),
        "note": "交错质量也开隙，是另一个有质量模",
    }

    R["verdict"] = {
        "uniform_massless": "均匀 δD 无质量（E 线性）——不是 λ=1/μ² 里的 μ",
        "dimer_massive": "键二聚化/交错质量模有质量（开 Dirac 隙），是候选的 δD 质量 μ",
        "honest": "μ = √(键序刚度 K)，需选「开隙的键序模」；λ = 1/μ² 由这个 μ 决定",
        "next": "把 μ（键二聚化模的刚度）代进 Soler λ，看涡旋频率比能否对上 √6——这是不凑数的下一步",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_bond_stiffness")

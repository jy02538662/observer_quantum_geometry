"""
共振判生死：对齐度/失配度 能否预测 耦合（电荷转移）？

共振机制（用户直觉）：「两个 D 的模流频率谱对齐 → 强耦合；失配 → 弱耦合」。

先澄清一个前提（否则「模流频率谱」没定义）：
  模流 σ_t = ρ^{it} · ρ^{-it} 定义在「观察者态 ρ」上，不是 D 上。
  尺度不变的 ρ = C/λ 与 D 无关 → 所有 D 的模流谱相同 → 共振退化（无信息）。
  要让共振非退化，必须让观察者态依赖 D。最自然选择 = 热态 ρ_D = e^{-βD}/Z，
  此时 log ρ_D = -βD + 常数 → 模流谱 = β × (D 谱)。即「模流共振」在热态观察者下
  退化为「D 能量谱重叠」。

所以本实验测的是：**两个 D 的能量谱「对齐/失配」能否预测它们的「耦合（电荷转移）」**。

判生死：
  - 若电荷转移被「全谱失配 Δ」预测 → 共振机制（对齐→耦合）有数值底座；
  - 若电荷转移只被「费米能级差 δE_F」预测（局部量），全谱 Δ 不预测 → 共振的
    「全谱对齐」版本判负（耦合由费米能级差驱动，不由全谱重叠驱动）。

Code: `py -m experiments.exp_resonance_alignment_coupling`
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# π 磁通 D（L×L）+ 局域缺陷
# ---------------------------------------------------------------------------
def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    return H

def pi_flux_defect(L, Vdef):
    H = pi_flux(L)
    H[0, 0] += Vdef
    return H

# ---------------------------------------------------------------------------
# 两个 D 的能量谱 + 三个「对齐/失配」量 + 费米能级差
# ---------------------------------------------------------------------------
def spectrum_measures(L, Vdef):
    N = L * L
    lam = np.sort(np.linalg.eigvalsh(pi_flux(L)))
    mu = np.sort(np.linalg.eigvalsh(pi_flux_defect(L, Vdef)))
    # ① 费米能级差（局部量，半满 = 第 N/2 个本征值）
    dEF = mu[N // 2 - 1] - lam[N // 2 - 1]
    # ② 全谱失配 Δ = Σ(λ-μ)²（全局量，= 1.29 的「能级失配」）
    Delta = float(np.sum((lam - mu) ** 2))
    # ③ 对齐度 A = 归一化内积（谱重叠，互补量）
    A = float(np.sum(lam * mu) / (np.linalg.norm(lam) * np.linalg.norm(mu)))
    return dEF, Delta, A

# ---------------------------------------------------------------------------
# 耦合 = 电荷转移（D-D 自反，复用 exp_doping_DD_selfadjoint 的机制）
# ---------------------------------------------------------------------------
def charge_transfer(L, Vdef, tp):
    N = L * L
    H1 = pi_flux(L)
    H2 = pi_flux_defect(L, Vdef)
    def idx(x, y):
        return (x % L) * L + (y % L)
    M = 2 * N
    Htot = np.zeros((M, M))
    Htot[:N, :N] = H1
    Htot[N:, N:] = H2
    for y in range(L):
        i = idx(L - 1, y); j = N + idx(0, y)
        Htot[i, j] -= tp; Htot[j, i] -= tp
    ev, evec = np.linalg.eigh(Htot)
    occ_idx = np.argsort(ev)[:M // 2]
    occ_vec = evec[:, occ_idx]
    C_left = occ_vec[:N, :] @ occ_vec[:N, :].conj().T
    C_right = occ_vec[N:, :] @ occ_vec[N:, :].conj().T
    nL = np.sum(np.linalg.eigvalsh(C_left))
    nR = np.sum(np.linalg.eigvalsh(C_right))
    return float(nL - nR)

# ---------------------------------------------------------------------------
# 扫描缺陷强度，看三个量各自与电荷转移的关系
# ---------------------------------------------------------------------------
L = 6
tp = 1.0
Vdefs = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0]

dEFs, Deltas, As, transfers = [], [], [], []
for Vdef in Vdefs:
    dEF, Delta, A = spectrum_measures(L, Vdef)
    transfer = charge_transfer(L, Vdef, tp)
    dEFs.append(round(dEF, 6)); Deltas.append(round(Delta, 6)); As.append(round(A, 8)); transfers.append(round(transfer, 6))

R["scan"] = {
    "缺陷强度 Vdef": Vdefs,
    "费米能级差 δE_F（局部量）": dEFs,
    "全谱失配 Δ=Σ(λ-μ)²（共振失配量）": Deltas,
    "对齐度 A=谱重叠（共振对齐量）": As,
    "电荷转移 n_transfer（耦合）": transfers,
}

# 相关性：Pearson 相关系数（粗糙判据）
def pearson(xs, ys):
    x = np.array(xs, dtype=float); y = np.array(ys, dtype=float)
    if np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])

R["correlation"] = {
    "Pearson(δE_F, n_transfer)": round(pearson(dEFs, transfers), 4),
    "Pearson(Δ, n_transfer)": round(pearson(Deltas, transfers), 4),
    "Pearson(对齐度A, n_transfer)": round(pearson(As, transfers), 4),
    "解读": "哪个相关系数 |r| 最大，就是电荷转移的驱动量（局部 δE_F vs 全局 Δ/A）",
}

# ---------------------------------------------------------------------------
# 判生死结论
# ---------------------------------------------------------------------------
R["verdict"] = {
    "测试了什么": "两个 D 的谱「对齐/失配」（共振量）能否预测「耦合」（电荷转移）",
    "前提澄清": "模流谱定义在观察者态 ρ 上；热态 ρ=e^{-βD} 下模流谱 = β×D 谱，故「模流共振」=「D 谱重叠」",
    "结果": "见 correlation——电荷转移是否被全谱失配 Δ 预测，还是只被费米能级差 δE_F 预测",
    "若 Δ 预测耦合": "共振「全谱对齐→耦合」有数值底座",
    "若只有 δE_F 预测耦合": "共振「全谱对齐→耦合」判负：耦合由费米能级差（局部量）驱动，不由全谱重叠（全局量）驱动——共振机制需收窄到「费米能级附近的谱」才可能成立",
}

report(R, "exp_resonance_alignment_coupling")

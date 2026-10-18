"""
「R 完满、D 欠缺」的公式化 + 程序验证（2026-10-15）

概念（承接 D-D 自反的深层解释）：
  - R 是完满的（连续本体，超有限 II₁ 因子，含所有视角）；
  - D = E(R) 是欠缺的（有限观察者的一次切割，只能切出有限视角）；
  - 「一个 D 总是欠缺 → 指向其他 D」是 E 有限性的结构推论，不加新公设。

要「落盘」，必须把它落成可算的公式，不能只凭感觉。本文给三个可算的「欠缺度」：

  ① 模流谱隙（欠缺 = 有隙）：
     观察者态 ρ = diag(λ₁…λ_N)，尺度不变 λ_i ∝ 1/i。
     模流频率 θ_ij = ln(λ_i/λ_j) = ln(j/i)，共 N² 个（有限、离散）。
     最小非零隙 Δ_min(N) = ln(N/(N-1)) > 0 —— 有隙 = 欠缺。
     Δ_min(N) ~ 1/N → 0 只在 N→∞（R 完满处）。

  ② 欠缺是结构必然（Wielandt–Wintner 版，符号）：
     有限 N×N 矩阵的模流生成元 [log ρ, ·] 恰有 N² 个本征值（有限），
     连续谱需要无限维 —— 单 D 永远够不到连续谱。

  ③ 观察切割的 Jones index（欠缺度 = 丢了多少视角）：
     条件期望 E : M_n⊗M_k → M_n（部分迹）丢因子 k²，
     Jones index [M_n⊗M_k : M_n] = k² → 相对 R（无限维）任何有限切割 index = ∞。

诚实边界：这三条都是「R 完满 D 欠缺」的可算版本，但「对齐度 A_μ → 耦合」这一环
（共振机制）仍没有结构等式，本脚本第 ④ 节只算 A_μ 本身、不算耦合（见诚实结论）。

Code: `py -m experiments.exp_R_complete_D_deficient`
"""
import numpy as np
from sympy import Matrix, Symbol, log as slog, simplify, symbols, Rational
from experiments._common import report

R = {}

# ===========================================================================
# ① 欠缺 = 有隙：单 D 的模流谱是离散的，最小隙 Δ_min(N) > 0
# ===========================================================================
def modular_frequencies(N):
    """ρ = diag(1/i)（尺度不变），模流频率 θ_ij = ln(λ_i/λ_j) = ln(j/i)。"""
    lam = 1.0 / np.arange(1, N + 1)          # λ_i = 1/i（未归一，不影响 θ）
    loglam = np.log(lam)
    theta = loglam[:, None] - loglam[None, :]  # θ_ij = ln λ_i - ln λ_j = ln(j/i)
    return theta

def min_gap(N):
    """最小非零 |θ_ij|（相邻对角元 = ln(N/(N-1))）。"""
    theta = modular_frequencies(N)
    theta_flat = np.abs(theta[np.abs(theta) > 1e-15])  # 去掉 θ_ii = 0
    return float(theta_flat.min())

Ns = [2, 3, 5, 10, 30, 100, 300, 1000]
gaps = [min_gap(N) for N in Ns]
theory_gap = [float(np.log(N / (N - 1))) for N in Ns]

R["1_deficient_gap"] = {
    "N": Ns,
    "最小模流谱隙 Δ_min(N)（数值）": [round(g, 10) for g in gaps],
    "理论 ln(N/(N-1))": [round(g, 10) for g in theory_gap],
    "数值=理论": [abs(a - b) < 1e-12 for a, b in zip(gaps, theory_gap)],
    "Δ_min(N) > 0 对所有有限 N（有隙=欠缺）": all(g > 0 for g in gaps),
    "Δ_min(N) ~ 1/N → 0（只在 N→∞ 完满）": f"gaps[N=1000] = {gaps[-1]:.6f} ≈ 1/1000",
}

# ===========================================================================
# ② 欠缺是结构必然：有限维模流生成元恰 N² 个本征值（符号）
# ===========================================================================
# 模流生成元 [log ρ, ·] 作用在矩阵单位 E_ij 上：
#   [log ρ, E_ij] = (log λ_i - log λ_j) E_ij = ln(λ_i/λ_j) E_ij
# 所以本征值 = {ln(λ_i/λ_j)}，恰 N² 个（有限、离散）。
n = 3
lam_sym = [Symbol(f"lambda_{i}") for i in range(n)]
# 符号验证：[log ρ, E_ij] 的对角展开
log_rho = Matrix.diag(*[slog(lam_sym[i]) for i in range(n)])
E = Matrix.zeros(n, n)
E[0, 1] = 1   # E_01
commutator = simplify(log_rho * E - E * log_rho)
expected = slog(lam_sym[0]) - slog(lam_sym[1])   # ln(λ_0/λ_1)
R["2_structural_symbolic"] = {
    "模流生成元 [log ρ, E_01]": str(commutator),
    "= (ln λ_0 - ln λ_1) E_01": str(simplify(commutator - expected * E)),
    "本征值 = {ln(λ_i/λ_j)}，恰 N² 个": n * n,
    "有限维 → 离散谱 → 有隙（结构必然，非数值巧合）": True,
}

# ===========================================================================
# ③ 观察切割的 Jones index：欠缺度 = 丢了多少视角
# ===========================================================================
# 部分迹 E : M_n⊗M_k → M_n（x0⊗x1 ↦ x0·τ(x1)），丢因子 k²。
# Jones index [M_n⊗M_k : M_n] = dim(M_n⊗M_k)/dim(M_n) = (nk)²/n² = k²。
n0 = 2
ks = [1, 2, 4, 8, 16, 64, 256]
indexes = [ (n0 * k) ** 2 // (n0 ** 2) for k in ks ]   # = k²
R["3_jones_index_deficiency"] = {
    "环境维度 k": ks,
    "Jones index [M_n⊗M_k : M_n] = k²": indexes,
    "欠缺度（丢的视角数）随 k 增长": indexes[-1],
    "相对 R（k→∞）任何有限切割 index=∞": "一个 D 永远欠缺（指向其他 D）",
    "assert_k_squared": indexes == [k * k for k in ks],
}

# ===========================================================================
# ④ 共振对齐度 A_μ（可算，但「→耦合」无结构等式，只算 A_μ 本身）
# ===========================================================================
def modular_density(N, n_bins=400, omega_max=10.0):
    """模流频率谱的（三角）密度：θ_ij = ln(j/i) 的直方图。"""
    theta = modular_frequencies(N)
    theta_flat = theta.flatten()
    theta_flat = theta_flat[~np.isclose(theta_flat, 0)]  # 去 θ_ii=0
    hist, edges = np.histogram(theta_flat, bins=n_bins, range=(-omega_max, omega_max), density=True)
    omega = 0.5 * (edges[:-1] + edges[1:])
    return omega, hist

def alignment(N1, N2, mu):
    """A_μ = ∫_0^μ ρ1(ω) ρ2(ω) dω（对齐度，重叠到尺度 μ）。"""
    w1, r1 = modular_density(N1)
    w2, r2 = modular_density(N2)
    # 插值到公共网格
    w = np.linspace(0, mu, 400)
    r1i = np.interp(w, w1, r1, left=0, right=0)
    r2i = np.interp(w, w2, r2, left=0, right=0)
    return float(np.trapezoid(r1i * r2i, w))

# 两个观察者（不同 N = 不同完满程度）的对齐度
align_same = alignment(100, 100, 10.0)
align_diff = alignment(100, 30, 10.0)
R["4_alignment_A_mu"] = {
    "对齐度 A_μ(ρ_N1, ρ_N2) 定义": "∫_0^μ ρ1(ω)ρ2(ω)dω（模流谱重叠）",
    "A_μ(100,100)（同完满度，自对齐）": round(align_same, 6),
    "A_μ(100,30)（不同完满度）": round(align_diff, 6),
    "A_μ 随观察者参数变（可算、可跑动）": align_same != align_diff,
    "诚实边界": "A_μ 可算，但「A_μ → 耦合 g」没有结构等式 g=g0·f(A_μ)，共振仍是类比不是推导",
}

# ===========================================================================
# 诚实结论
# ===========================================================================
R["honest_conclusion"] = {
    "R 完满 D 欠缺的可算版本": {
        "① 欠缺=有隙": "单 D 模流谱最小隙 Δ_min(N)=ln(N/(N-1))>0，~1/N",
        "② 结构必然": "有限维模流生成元恰 N² 个本征值（符号），连续谱需无限维",
        "③ Jones index": "切割 E 丢 k² 视角，相对 R 任何有限切割 index=∞",
    },
    "完满在极限处": "Δ_min(N)→0 只在 N→∞（II₁ 非迹权给连续谱，1.3 §2 已证）",
    "与已有结果的关系": "「欠缺=有隙/完满=无隙」是 1.3 §1.3「谱隙↔关联长程」的重述，新增的是「欠缺是结构必然 + 一个 D 指向其他 D」的概念澄清",
    "还没算的": "「对齐度 A_μ → 耦合」这一步（共振机制的核心）仍无结构等式，本脚本不假装算出了它",
}

report(R, "exp_R_complete_D_deficient")

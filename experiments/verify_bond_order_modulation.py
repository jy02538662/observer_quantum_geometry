"""
验证「内生无序」路径三：确定性键序反馈 D = D_0 + λρ[D] 是否自发产生
空间调制（CDW/对称破缺）？

背景（2026-10-16）：
  「内生无序」三条候选路径里，路径一（多 D 集体效应）= 手放推到上一层；
  路径二（模流涨落）= 被幺正性杀掉（模流 ρ^{it}·ρ^{-it} 是幺正的，谱不变，
  无 Lyapunov 增长，无「湍流」）。路径三是唯一有真实物理的：
  确定性反馈 D = D_0 + λρ[D]（自洽 Hartree-Fock 型）能否自发破缺平移对称。

  关键物理：半填充 π 磁通 = Dirac 半金属。Dirac 半金属对弱相互作用是
  「半金属稳定」（DOS 在 Dirac 点消失，需有限阈值 U_c 才能开隙）。
  所以问题 = 反馈在 λ 多大时会自发产生 CDW/调制。

  本脚本：N×N π 磁通环面，自洽迭代 D <- D_0 + λρ[D]（ρ = 价带投影），
  加微小随机种子破对称，测「调制振幅」（ρ_ii 偏离均匀 1/2 的程度）随 λ 变化。

Code: `py -m experiments.verify_bond_order_modulation`
"""
import numpy as np

print("=" * 72)
print("验证：键序反馈 D = D_0 + λρ[D] 是否自发产生调制（CDW）？")
print("=" * 72)

N = 8
n = N * N  # 64 sites
nocc = n // 2  # half filling = 32


def idx(x, y):
    return (x % N) * N + (y % N)


# π-flux H_0（t=1）
H0 = np.zeros((n, n))
for x in range(N):
    for y in range(N):
        i = idx(x, y)
        j = idx(x + 1, y)
        H0[i, j] -= 1.0
        H0[j, i] -= 1.0
        j = idx(x, y + 1)
        phase = (-1.0) ** x
        H0[i, j] -= phase
        H0[j, i] -= phase


def gs_rho(D):
    """基态密度矩阵 ρ = 价带投影（最低 nocc 个本征态）。"""
    evals, evecs = np.linalg.eigh(D)
    occ = evecs[:, :nocc]
    return occ @ occ.T, evals


def iterate(lam, niter=200, seed_amp=1e-3):
    """自洽迭代 D <- D_0 + λρ，从 H0 + 微小随机扰动开始。"""
    rng = np.random.default_rng(42)
    # 微扰破对称（对角随机）
    D = H0.copy()
    for i in range(n):
        D[i, i] += seed_amp * rng.standard_normal()
    for _ in range(niter):
        rho, evals = gs_rho(D)
        Dnew = H0 + lam * rho
        D = Dnew
    rho, evals = gs_rho(D)
    return rho, evals


print("\n[调制振幅 vs λ]（modulation = max|ρ_ii - 1/2|，0 = 均匀无破缺）")
print(f"{'λ':>6} {'调制振幅':>10} {'能隙':>10} {'说明':>24}")

for lam in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
    rho, evals = iterate(lam)
    diag = np.diag(rho)
    mod = np.max(np.abs(diag - 0.5))
    # 能隙 = E[nocc] - E[nocc-1]（半填充处）
    gap = evals[nocc] - evals[nocc - 1]
    if mod < 1e-6:
        note = "均匀（无破缺）"
    elif mod < 0.05:
        note = "弱调制（阈值附近）"
    else:
        note = "强调制（CDW 破缺）"
    print(f"{lam:>6.1f} {mod:>10.5f} {gap:>10.5f} {note:>24}")

print("\n" + "=" * 72)
print("结论")
print("=" * 72)
print("若 λ 超过某个有限阈值才出现调制（Dirac 半金属对弱反馈稳定），")
print("则反馈「能」内生对称破缺，但：")
print("  1. 需要有限阈值 λ_c（不是任意小 λ）；")
print("  2. 产生的是「有序」调制（CDW，周期性），不是「随机」无序；")
print("  3. 周期调制开的是「能隙」（能带绝缘体/CDW），不是 Anderson 局域化；")
print("  4. 这是自洽平均场（单粒子），不是真正的多体 Hubbard U（Mott）。")

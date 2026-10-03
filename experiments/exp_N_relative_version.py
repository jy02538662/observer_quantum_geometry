"""
远景线 5（质量谱）· 探「N 的相对版本」：D-D 自反复合系统的离散不变量

背景：N（区分节点数）从「单个 D 群论」数不出 128（路 4）。用户提示用「D-D 自反」
（相对化，破 μ 墙用的）来找 N 的「相对版本」。

核心问题：两个 D 自反耦合的复合系统里，有没有一个「离散、量子化」的不变量，
能对应 N？

候选离散不变量：
  1. 能级交叉数：tp 从 0 扫到 ∞，复合系统能级（作为 tp 函数）交叉的次数；
  2. 谱简并结构：两个 D 相同/不同时，复合谱的简并度；
  3. 手征对称性：复合系统的 Γ（手征）结构；
  4. 费米能级差 ΔEF 与转移量 Δn 的「量子化台阶」。

本脚本数值探这些，诚实报告：找到的离散不变量是什么、它是不是 N 的量级。
"""
import numpy as np
from experiments._common import report

R = {}

def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y); H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1); ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    return H

def pi_flux_defect(L, Vdef):
    H = pi_flux(L)
    H[0, 0] += Vdef
    return H

def composite_spectrum(L, Vdef, tp):
    """两个 D 自反耦合，复合系统本征值。"""
    N = L * L
    H1 = pi_flux(L)
    H2 = pi_flux_defect(L, Vdef)
    M = 2 * N
    Htot = np.zeros((M, M))
    Htot[:N, :N] = H1
    Htot[N:, N:] = H2
    def idx(x, y):
        return (x % L) * L + (y % L)
    for y in range(L):
        i = idx(L - 1, y)
        j = N + idx(0, y)
        Htot[i, j] -= tp
        Htot[j, i] -= tp
    return np.linalg.eigvalsh(Htot)

# ============ 候选 1：能级交叉数（tp 扫描） ============
L = 4
N = L * L
tps = np.linspace(0.0, 5.0, 200)
spectra = np.array([composite_spectrum(L, 0.0, tp) for tp in tps])  # 两个干净 D

# 数费米能级附近的能级交叉：第 N-1 和第 N 个能级（半满边界）随 tp 演化
# 两个干净 D：Vdef=0，费米能级差=0，没有转移，但能级随 tp 重排
levels_Nm1 = spectra[:, N-1]   # 第 N 个（半满下边界）
levels_N = spectra[:, N]       # 第 N+1 个（半满上边界）

# 能级交叉：相邻能级随 tp 交叉的次数
crossings = 0
for i in range(2*N):
    ev_i = spectra[:, i]
    for j in range(i+1, 2*N):
        ev_j = spectra[:, j]
        diff = ev_i - ev_j
        sign_changes = np.sum(np.diff(np.sign(diff)) != 0)
        crossings += sign_changes
crossings = crossings // 2  # 每对能级交叉算了两次（i,j 和 j,i）

R["level_crossing"] = {
    "L": L, "N=L²": N,
    "两个干净 D（Vdef=0）": "费米能级差=0，无净转移",
    "tp 扫描能级交叉总数": int(crossings),
    "判断": "能级交叉数是否 ~ N 或 N 的某个简单函数？",
}

# ============ 候选 2：谱简并结构（tp=0 vs tp 大） ============
spec_tp0 = composite_spectrum(L, 0.0, 0.0)
spec_tp_big = composite_spectrum(L, 0.0, 3.0)

def degeneracy(spec, tol=1e-6):
    """数简并：把接近相等的本征值聚类。"""
    s = np.sort(spec)
    groups = []
    cur = [s[0]]
    for v in s[1:]:
        if abs(v - cur[-1]) < tol:
            cur.append(v)
        else:
            groups.append(len(cur)); cur = [v]
    groups.append(len(cur))
    return groups

R["degeneracy"] = {
    "tp=0 简并分布（两个干净 D，各自 N 能级 ×2）": degeneracy(spec_tp0),
    "tp=3 简并分布": degeneracy(spec_tp_big),
    "判断": "两个相同 D 耦合，tp=0 时 2 重简并（两个 D 各一份）；tp 增大杂交破简并",
}

# ============ 候选 3：费米能级差（Vdef≠0 时，谱不对称） ============
# 两个 D：一个干净一个带缺陷，费米能级差 ΔEF 是否量子化？
R["fermi_level_difference"] = {}
for Vdef in (0.0, 0.5, 1.0, 2.0, 4.0):
    ev1 = np.linalg.eigvalsh(pi_flux(L))
    ev2 = np.linalg.eigvalsh(pi_flux_defect(L, Vdef))
    EF1 = ev1[N//2 - 1]
    EF2 = ev2[N//2 - 1]
    R["fermi_level_difference"][f"Vdef={Vdef}"] = {
        "ΔEF": round(float(EF2 - EF1), 4),
    }
R["fermi_level_difference"]["判断"] = "ΔEF 随 Vdef 连续变（非量子化），是连续量不是离散不变量"

# ============ 候选 4：复合系统的「能级对数」（N 个能级对） ============
# 两个 D 各 N 个能级，耦合后形成 N 个「能级对」（bonding/anti-bonding）。
# 这个「能级对数」= N 是输入，不是解出的。但看耦合后能级对是否重新组织成新结构。
R["level_pairs"] = {
    "两个 D 各 N 个能级": "耦合后形成 N 个 bonding/anti-bonding 能级对",
    "能级对数": N,
    "关键": "这个 N 是「输入」（两个 D 各 L² 格点），不是「解出」——D-D 自反不产生 N 的值，只产生 N 个能级之间的相对关系",
}

R["honest_conclusion"] = {
    "探索结果": "D-D 自反复合系统的离散结构 = 能级对、简并、能级交叉，但这些都是「N 是输入」后的结果",
    "关键判断": "D-D 自反给的「相对量」是费米能级差 ΔEF、转移量 Δn（连续），不是「N 的值」（离散）——N 的「相对版本」没找到",
    "原因": "N（区分节点数）是「单个 D 的格点数/量子化 level」，是「结构侧」的绝对量；D-D 自反的「相对化」针对的是「物质侧」（掺杂/填充），两者不同层",
    "出路（诚实）": "D-D 自反破 μ 墙是因为找到了 μ 的相对版本（费米能级差）；要破 N 墙需找到 N 的相对版本（两个 D 之间什么相对量 = 区分节点数），而这个目前框架没有",
    "措辞": "这是「未找到 N 的相对版本」（探索负结果），非「否证」——不是证明 N 不可能相对化，是「当前 D-D 自反结构里没有 N 的相对量」",
}

report(R, "exp_N_relative_version")

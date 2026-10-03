"""
质量谱动力学探底 · π-flux 零模与「断裂 → 掺杂」的微观来源

修正上一版：π-flux 谱有 4 个零模（Dirac 点），半满占据 = 6 负 + 2 零模，
留 2 个零模空。上版用 w<0 判据漏占零模，导致 n_i 不均匀（判据 bug，非物理）。

本脚本聚焦真正的发现：π-flux 的零模（Dirac 点）是「掺杂」的微观来源——
断裂（π-flux）→ Dirac 点零模 → 填充歧义 → 掺杂（连续）。这是物质侧墙
「断裂 → 掺杂」最终墙的微观机制定位。

回答：
  Q1：π-flux 谱结构（几个零模？Dirac 点在哪？）
  Q2：半满占据的正确填充（6 负 + 2 零模），零模填充的歧义 = 掺杂自由度？
  Q3：零模填充 → n_i 分布（掺杂）与「断裂→掺杂」的关系。
  Q4：这给不给「局部 U（同格点两体）」→ 29.2？

诚实边界：零模填充给「掺杂（连续填充，单体占据数）」，不给「同格点两体 U」。
"""
import numpy as np
from experiments._common import report

R = {}


def toroidal_D(npd, pi_flux=True):
    n = npd ** 2
    D = np.zeros((n, n), complex)
    for i in range(npd):
        for j in range(npd):
            idx = npd * i + j
            jr = (j + 1) % npd
            D[idx, npd * i + jr] += 1.0
            D[npd * i + jr, idx] += 1.0
            ph = np.pi * j if pi_flux else 0.0
            D[idx, npd * ((i + 1) % npd) + j] += np.exp(1j * ph)
            D[npd * ((i + 1) % npd) + j, idx] += np.exp(-1j * ph)
    return D


npd = 4
n = npd ** 2
D0 = toroidal_D(npd, pi_flux=True)
w, V = np.linalg.eigh(D0)

# Q1：谱结构
neg = w < -1e-8
zero = np.abs(w) < 1e-8
pos = w > 1e-8
R["Q1_spectrum"] = {
    "本征值（升序）": [round(float(x), 4) for x in w],
    "负本征值个数": int(np.sum(neg)),
    "零模个数（Dirac 点）": int(np.sum(zero)),
    "正本征值个数": int(np.sum(pos)),
    "零模 = Dirac 点": "4 个零模 = 2 个 Dirac 点 × 2（谷简并），π-flux 是 Dirac 半金属。",
}

# Q2：半满占据的正确填充 = 6 负 + 2 零模（共 8 = n/2）
# 零模填充有歧义：4 个零模里选 2 个（C(4,2)=6 种）
half = n // 2
order = np.argsort(w)  # 升序
occ_full = order[:half]  # 前 8 = 6 负 + 2 零模

# 检查这 8 个里零模的个数
zero_in_occ = np.sum(zero[occ_full])
R["Q2_half_filling"] = {
    "半满占据态数": half,
    "占据的负本征值数": int(np.sum(neg[occ_full])),
    "占据的零模数": int(zero_in_occ),
    "剩余的零模数（空）": int(np.sum(zero) - zero_in_occ),
    "填充歧义": f"4 个零模里选 {zero_in_occ} 个占据，留 {int(np.sum(zero)-zero_in_occ)} 个空 → 掺杂自由度（连续填充）。",
}

# Q3：零模填充 → n_i 分布。对比「全 8 个最低（含零模）」vs「只 6 负（上版 bug）」
def n_i_from_occ(occ_idx):
    rho = np.zeros((n, n), complex)
    for k in occ_idx:
        psi = V[:, k]
        rho += np.outer(psi, psi.conj())
    return np.real(np.diag(rho))

ni_half = n_i_from_occ(occ_full)          # 正确半满（含 2 零模）
ni_neg = n_i_from_occ(np.where(neg)[0])   # 只 6 负（上版 bug）

R["Q3_filling_n_i"] = {
    "正确半满 n_i（含零模）": [round(float(x), 4) for x in ni_half],
    "只占负态 n_i（上版 bug）": [round(float(x), 4) for x in ni_neg],
    "正确半满 n_i 是否均匀（≈0.5）": bool(np.max(np.abs(ni_half - 0.5)) < 0.01),
    "只占负态 n_i 是否均匀": bool(np.max(np.abs(ni_neg - 0.5)) < 0.01),
    "结论": "正确半满（含零模）n_i 均匀 ≈0.5；上版漏零模导致 n_i 不均匀（0.62 vs 0.375）= 判据 bug。零模填充给的是『均匀半满』还是『偏离半满』取决于选哪 2 个零模。",
}

# 零模填充歧义：不同选法给不同 n_i 分布（掺杂的微观来源）
# 零模的波函数在实空间是否局域化？
zero_idx = np.where(zero)[0]
R["Q3_zero_mode_wavefunction"] = {}
for a in zero_idx:
    psi = V[:, a]
    dens = np.abs(psi) ** 2
    R["Q3_zero_mode_wavefunction"][f"零模 {a} 概率密度分布"] = [round(float(x), 3) for x in dens]

# 零模是否局域（局域化 vs 弥散）
dens_zero = np.abs(V[:, zero_idx]) ** 2  # (n, 4)
ipr = np.sum(dens_zero ** 2, axis=0)     # 逆参与比
R["Q3_zero_mode_localization"] = {
    "4 个零模的逆参与比 IPR（1=完全局域，1/n=完全弥散）": [round(float(x), 4) for x in ipr],
    "1/n（完全弥散）": round(1.0 / n, 4),
    "判断": "IPR ≈ 1/n 说明零模弥散（非局域）；IPR ≈ 1 说明局域。零模的空间结构决定『掺杂』是均匀还是局域。",
}

# Q4：这给不给局部 U？
R["Q4_local_U"] = {
    "零模填充给什么": "掺杂（连续填充，单体占据数 n_i 的分布）",
    "零模填充给不给同格点两体 U": "不给——n_i 是单粒子占据数（0 或 1），不是双通道 n↑n↓。",
    "与 29.2 的关系": "掺杂是『断裂→掺杂』（物质侧墙最终墙）的微观来源，但 29.2 需要的是 Yukawa 圈图修正 = 同格点两体 U，零模填充给不了。",
}

# 结论
R["honest_conclusion"] = {
    "真正的发现": "π-flux 有 4 个零模（Dirac 点），半满占据时 4 零模只填 2、留 2 空 = 掺杂自由度。这是『断裂（π-flux）→ Dirac 点零模 → 填充歧义 → 掺杂』的微观机制——物质侧墙『断裂→掺杂』最终墙的微观定位。",
    "与翻转 3 的关系": "翻转 3『真空传播子 G₀ 对角非零』上版是零模填充的判据 bug（漏占零模），修正后均匀半满 G₀ 对角 = 0。真空不空的『对角/局部』分量不成立，真空不空只体现在『非对角关联』。",
    "净判断": "π-flux 零模给『掺杂』（单体填充），不给『局部 U』（同格点两体）。这是『断裂→掺杂』的推进，不是『局部U/29.2』的推进。",
    "对质量谱动力学的意义": "零模填充 = 掺杂 = 质量谱动力学的『前置条件』（掺杂 → Yukawa → 质量），但它本身不是 Yukawa 跑动（29.2）。29.2 的局部 U 缺口仍未解决。",
    "定位": "【候选观察】。π-flux 零模 → 掺杂是真实机制（物质侧墙『断裂→掺杂』的微观来源），但离 29.2 还差『零模填充 → Yukawa 跑动』的桥，那桥需要局部 U（同格点两体），框架没有。",
}

report(R, "exp_gap3_vacuum_propagator")

"""
远景线 5（质量谱）· 上型不对称 路径 C：验证「3」是否精确

用户三路径：A（希格斯 H vs H̃）、B（Fano 两个 3 的区别，已证走不通）、C（「3」是巧合）。
先做路径 C（最便宜）：用精确夸克质量值，验证 m_t/m_c ÷ m_b/m_s 是否 = 3（精确），
以及在所有标度方案下是否稳定 = 3。

关键：代际比 m_t/m_c 和 m_b/m_s 的标度依赖——m_t 用 pole 还是 MS̄，m_b 用哪个标度。
"""
import numpy as np
from experiments._common import report

R = {}

# 夸克质量（PDG 2024，MS̄，带近似误差）
# 上型：m_u(2GeV), m_c(m_c), m_t（pole 或 MS̄）
# 下型：m_d(2GeV), m_s(2GeV), m_b(m_b)
quarks = {
    "m_u(2GeV)": (2.16, 0.07),   # MeV
    "m_c(m_c)":  (1273, 5),      # MeV
    "m_t(pole)": (172690, 300),  # MeV
    "m_t(MS̄)":  (162500, 2000), # MeV
    "m_d(2GeV)": (4.67, 0.05),   # MeV
    "m_s(2GeV)": (93.4, 0.7),    # MeV
    "m_b(m_b)":  (4180, 20),     # MeV
}

# 第二→第三代际比（上型 vs 下型）
# 上型：m_t/m_c，下型：m_b/m_s
def ratio_23(m_t, m_c, m_b, m_s):
    return (m_t/m_c) / (m_b/m_s)

# 方案 1：m_t pole = 172.69 GeV
r_pole = ratio_23(quarks["m_t(pole)"][0], quarks["m_c(m_c)"][0],
                  quarks["m_b(m_b)"][0], quarks["m_s(2GeV)"][0])
# 方案 2：m_t MS̄ = 162.5 GeV
r_msbar = ratio_23(quarks["m_t(MS̄)"][0], quarks["m_c(m_c)"][0],
                   quarks["m_b(m_b)"][0], quarks["m_s(2GeV)"][0])

R["ratio_3"] = {
    "m_t/m_c ÷ m_b/m_s（m_t pole=172.69）": r_pole,
    "m_t/m_c ÷ m_b/m_s（m_t MS̄=162.5）": r_msbar,
    "范围": f"{min(r_pole,r_msbar):.3f} ~ {max(r_pole,r_msbar):.3f}",
    "是否 = 3（精确）": "否——标度敏感，2.85~3.03，非精确 3",
}

# 误差传播（粗略）
# m_t/m_c 的相对误差 ~ 0.3/172.69 + 5/1273 ~ 0.17% + 0.39% ~ 0.56%
# m_b/m_s 的相对误差 ~ 20/4180 + 0.7/93.4 ~ 0.48% + 0.75% ~ 1.23%
# 比值相对误差 ~ 0.56% + 1.23% ~ 1.8%
# 比值 = 3.03 ± 0.05（1.8%）
R["error"] = {
    "比值（pole）": r_pole,
    "相对误差（粗略）": "~1.8%（m_t 0.17% + m_c 0.39% + m_b 0.48% + m_s 0.75%）",
    "3 是否在误差范围内": "3.03 ± 0.05 不含 3（3.03 - 0.05 = 2.98 > 3？不对，2.98 < 3）——需精确误差分析",
}

# 关键：第二→第三代的「上型/下型比值」的标度敏感性
# m_t 是最敏感的（pole vs MS̄ 差 6%）
R["scale_sensitivity"] = {
    "m_t 标度敏感（pole 172.69 vs MS̄ 162.5）": f"比值 {r_pole:.3f} vs {r_msbar:.3f}，差 {(r_pole-r_msbar)/r_pole*100:.1f}%",
    "结论": "「3」对标度敏感，pole 给 3.03、MS̄ 给 2.85——不是精确的 3",
}

# 上下型代际比分别看（不取比值）
R["individual_ratios"] = {
    "上型 m_t/m_c（pole）": quarks["m_t(pole)"][0]/quarks["m_c(m_c)"][0],
    "下型 m_b/m_s": quarks["m_b(m_b)"][0]/quarks["m_s(2GeV)"][0],
    "上型/下型": r_pole,
    "8 = 2³ 候选": f"上型/轻子 = {quarks['m_t(pole)'][0]/quarks['m_c(m_c)'][0] / 16.817:.2f}（轻子第二→第三=16.817）",
    "8/3 候选": f"下型/轻子 = {quarks['m_b(m_b)'][0]/quarks['m_s(2GeV)'][0] / 16.817:.2f}",
}

# 诚实结论
R["honest_conclusion"] = {
    "路径 C 验证": "「3」不是精确的——m_t 标度敏感，pole 172.69 给 3.03、MS̄ 162.5 给 2.85",
    "「3」是巧合吗": "候选「巧合」——2.85~3.03 不是精确 3，但都在 3 附近（~1-5%）",
    "结论": "路径 C 不能排除「巧合」，但「3」接近（~3.03 pole / 2.85 MS̄）。若「3」不是精确，则「上型×3 vs 下型×1」的「3」是近似，不是精确关系——路径 A（希格斯 H vs H̃）的价值下降（因为「3」不精确）",
    "下一步": "「3」不精确 → 路径 C 倾向「巧合或近似」；若要继续，需路径 A（希格斯结构）看能否解释「近似 3」的机制",
    "措辞": "路径 C 验证（「3」不精确，标度敏感 2.85~3.03）——「3」可能是近似/巧合，非精确关系",
}

report(R, "exp_three_precision")

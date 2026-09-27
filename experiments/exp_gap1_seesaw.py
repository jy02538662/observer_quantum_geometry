"""
v11 缺口 1 · seesaw 的 1.6 倍因子：精确公式（1/2 因子）+ Y_ν 标定

问题：框架 m_ν = v²/M_R 给 81 meV，观测最重中微子 ~50 meV，差 1.6 倍。
精确化：
  type-I seesaw 精确公式 m_ν = m_D²/M_R = (Y_ν v/√2)²/M_R = Y_ν² v²/(2M_R)
  （代码 m_ν = v²/M_R 丢了 1/2 和 Y_ν²）

链条：
  M_R = M_P/N_int²（1生2：N_int 幂次=区分次数，M_R∝m²）
  m_ν = Y_ν² v²/(2 M_R) = Y_ν² v² N_int²/(2 M_P)
"""
import numpy as np
from experiments._common import report

R = {}

v = 246.22          # GeV 电弱 VEV
M_P = 1.22e19       # GeV
N_int = 128         # 2^7（三个 Z₂ → Fano → 7 → 128）

M_R = M_P / N_int**2

# 代码旧公式（无 1/2）
m_nu_old = v**2 / M_R          # = 81.4 meV
# 精确 seesaw 公式（补 1/2）
m_nu_exact = v**2 / (2 * M_R)  # = 40.7 meV

R["M_R"] = {
    "M_R = M_P/N_int²": f"{M_R:.2e} GeV = {M_R/1e15:.2f}×10^15 GeV（GUT 标度）",
}

R["two_formulas"] = {
    "代码旧公式 m_ν = v²/M_R（无 1/2）": f"{m_nu_old*1e12:.1f} meV",
    "精确 seesaw m_ν = v²/(2M_R)（补 1/2）": f"{m_nu_exact*1e12:.1f} meV",
    "差": f"{m_nu_old/m_nu_exact:.1f} 倍（=2，即缺的 1/2）",
}

# 观测约束
dmsq21 = 7.5e-5   # eV² 太阳
dmsq31 = 2.5e-3   # eV² 大气
m3_min = np.sqrt(dmsq31) * 1e3   # 正常序最重 m3 下限（m1→0），meV
R["obs"] = {
    "Δm²_31 = 2.5e-3 eV²": f"最重中微子 m3 下限 = √(Δm²_31) = {m3_min:.1f} meV（正常序，m1→0）",
    "宇宙学上限": "Σm_ν < 120 meV",
    "观测范围": "最重中微子 ~ 50 meV（下限）",
}

# 匹配观测最重 m3=50 meV，反推 Y_ν
m_nu_target = 50e-3 * 1e-9  # 50 meV = 5e-11 GeV
Y_nu2 = m_nu_target / m_nu_exact
Y_nu = np.sqrt(Y_nu2)
R["Y_nu"] = {
    "Y_ν² = m3/(v²/(2M_R))": f"{Y_nu2:.3f}",
    "Y_ν": f"{Y_nu:.3f}",
    "1.6 因子 = 2/Y_ν²": f"{2/Y_nu2:.2f}（代码 81 meV / 观测 50 meV = 81/50 = 1.62）",
}

R["honest_conclusion"] = {
    "精确化": "代码 m_ν=v²/M_R 丢了 type-I seesaw 的 1/2 因子。补上后 m_ν = Y_ν² v²/(2M_R) = 40.7 meV × Y_ν²。",
    "1.6 的来源": "1.6 = 2/Y_ν² =（缺 1/2 的代码 bug）×（Y_ν≈1.1 的标定）。81 meV（无 1/2）/50 meV（观测）=1.62。",
    "Y_ν 是标定": "Y_ν ≈ 1.1 是 O(1) 标定输入——「无绝对尺度」⟹ 绝对耦合（Yukawa 系数）是标定，不是框架能推的纯数。与轻子质量谱（1.14）的「无自由参数」不同，seesaw 的 Y_ν 是标定。",
    "1生2 不受影响": "M_R = M_P/N_int²（两次区分，M_R∝m²）的幂次结构不受 1/2 系数影响（1/2 是 seesaw 公式系数，非幂次）。",
    "措辞": "三档「未证实」——Y_ν 是标定输入（默认档），非「否证」（没证明 Y_ν 绝不可能有框架理由），也非「重言式」。1/2 因子是代码精确化（补 bug），非结构推进。",
}

report(R, "exp_gap1_seesaw")

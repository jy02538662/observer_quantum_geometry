"""
远景线 5（质量谱）· 精确化：质量 = 频率比（操作定义 ↔ 模流）

用户洞察（关键）：质量测量的「操作定义」本质上是「频率比」——
  彭宁阱：ω_c = qB/m（回旋频率）vs 拉莫尔频率；
  μ子偶素：超精细跃迁频率；
  τ：衰变运动学边缘（赝质量）。
而框架的「频率」= 模流（时间演化）σ_t = ρ^{it}，频率 = log ρ 的谱。

一元论精确化：
  (1) 质量 ∝ 能量频率（E = mc² = ℏω ⟹ m ∝ ω）；
  (2) 质量比 = 频率比（无量纲，框架无绝对尺度只能预言比值——与 CODATA 用比值定义常数自洽）；
  (3) 「代」= S₃ 共轭类（阶 1,2,3）→ 阶是「离散模式」，映射到「连续频率」（模流）= 离散→连续（一元论翻转）。

本脚本验证「阶 → 频率」的映射结构，诚实标注「纯指数对不上（需修正）」。
"""
import numpy as np
from experiments._common import report

R = {}

# 观测数据（CODATA 2022，比值）
m_ratio_mu_e = 206.7682830      # m_μ/m_e
m_ratio_tau_e = 3477.0          # m_τ/m_e（约）

# 对数（质量 ∝ 频率 ⟹ 质量比的对数 = 频率比的对数）
ln_mu_e = np.log(m_ratio_mu_e)
ln_tau_e = np.log(m_ratio_tau_e)

R["observed"] = {
    "m_μ/m_e": m_ratio_mu_e,
    "m_τ/m_e": m_ratio_tau_e,
    "ln(m_μ/m_e)": float(ln_mu_e),
    "ln(m_τ/m_e)": float(ln_tau_e),
    "ln(m_τ/m_e)/ln(m_μ/m_e)": float(ln_tau_e / ln_mu_e),
}

# 纯指数假设 m_n ∝ e^{c·n}（n=阶 1,2,3）：
#   ln(m_μ/m_e) = c（阶 2-1=1），ln(m_τ/m_e) = 2c（阶 3-1=2）
#   预测 ln(m_τ/m_e)/ln(m_μ/m_e) = 2
R["pure_exponential_test"] = {
    "假设": "m_n ∝ e^{c·n}（阶 n）",
    "预测 ln(m_τ/m_e)/ln(m_μ/m_e)": 2.0,
    "实际": float(ln_tau_e / ln_mu_e),
    "结论": "实际 1.53 ≠ 2 —— 纯指数不对（m_τ/m_e ≠ (m_μ/m_e)²），需 S₃ 更深结构的修正",
}

# 幂律假设 m_n ∝ n^k：
#   ln(m_μ/m_e) = k·ln(2)，ln(m_τ/m_e) = k·ln(3)
#   k = ln(m_μ/m_e)/ln(2)，预测 m_τ/m_e = 3^k
k = ln_mu_e / np.log(2)
pred_tau_power = 3**k
R["power_law_test"] = {
    "假设": "m_n ∝ n^k（阶 n）",
    "k = ln(m_μ/m_e)/ln2": float(k),
    "预测 m_τ/m_e = 3^k": float(pred_tau_power),
    "实际 m_τ/m_e": m_ratio_tau_e,
    "结论": f"预测 {pred_tau_power:.0f} vs 实际 3477 —— 纯幂律也不对（差 {(pred_tau_power/m_ratio_tau_e):.1f} 倍）",
}

R["honest_conclusion"] = {
    "一元论方向": "质量 = 频率比（操作定义）+ 频率 = 模流（时间），「阶 → 频率」= 离散→连续（一元论翻转）——方向是对的，精确化了「质量比 = 模流频率比」",
    "诚实发现": "纯指数（m_n∝e^{c·n}）和纯幂律（m_n∝n^k）都对不上 206.77 和 3477 —— 说明「阶 → 频率」的映射不是「单一阶的不变量」，需要「S₃ 的更深结构」（中心化子阶 6,2,3 / 类大小 1,3,2）的组合",
    "关键信号": "ln(m_τ/m_e)/ln(m_μ/m_e) = 1.53（不是纯指数的 2）——这个「1.53」的偏离，可能正是 S₃ 共轭类「非均匀结构」的指纹",
    "下一步": "「阶 → 频率」需要「两个不变量」的组合（阶 + 中心化子/类大小），而不是「单变量函数」——这是「猜数字」和「结构逼出」的分界",
    "措辞": "精确化方向（质量=频率比），非「已推出数值」——具体数值仍需「阶+中心化子」的二元结构",
}

report(R, "exp_mass_frequency")

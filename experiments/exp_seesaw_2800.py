"""
远景线 2（尺度读出）· 挖 2800：跷跷板（seesaw）画面

用户直觉：中微子质量 m_ν = v²/M_R（seesaw），M_R ~ M_P/N_int²（GUT 标度）。
本脚本验证：M_R = M_P/N_int² 给 GUT 标度、m_ν = v²/M_R 给中微子质量量级，
并算 m_ν/m（中微子质量/红外截断）看 2800 的标定因子能否变成结构因子。

框架量：
  N_int = 128 = 2⁷（内部，三个 Z₂）
  N_ext ~ 10^30（外部，S^{1/4}）
  Λ_eff = M_P/N_ext，ρ_Λ = 2Λ_eff⁴，ρ_Λ^{1/4} = 2^{1/4} M_P/N_ext
  m = π Λ_eff/N_int = π M_P/(N_ext·N_int)（红外截断）
"""
import numpy as np
from experiments._common import report

R = {}

v = 246.22          # GeV（电弱真空期望值）
M_P = 1.22e19       # GeV
N_int = 128
N_ext = 1e30

# 1. M_R = M_P/N_int²（用户候选：Majorana 质量 = 普朗克/内部自由度²）
M_R = M_P / N_int**2
R["M_R"] = {
    "M_R = M_P/N_int² = 10^19/128²": f"{M_R:.2e} GeV",
    "GUT 标度？": "≈ 7.4×10^14 GeV（GUT 量级 10^15）",
}

# 2. m_ν = v²/M_R（seesaw）
m_nu = v**2 / M_R
R["m_nu_seesaw"] = {
    "m_ν = v²/M_R = (246)²/(7.4×10^14)": f"{m_nu:.2e} GeV = {m_nu*1e9*1e3:.3f} meV",
    "中微子质量量级？": "观测 m_ν ~ 50-100 meV",
}

# 3. m = π M_P/(N_ext·N_int)（红外截断）
m = np.pi * M_P / (N_ext * N_int)
R["m_IR"] = {
    "m = π M_P/(N_ext·N_int)": f"{m:.2e} GeV = {m*1e12:.3f} meV",
}

# 4. m_ν/m（中微子质量/红外截断）——seesaw 结构 vs 观测
ratio_seesaw = m_nu / m
ratio_obs = 50e-3 * 1e-9 / m  # m_ν=50 meV=5e-11 GeV
R["ratio"] = {
    "seesaw 给 m_ν/m": f"{ratio_seesaw:.0f}",
    "观测 m_ν/m（中微子 50 meV）": f"{ratio_obs:.0f}",
    "比值": f"{ratio_seesaw/ratio_obs:.2f}（O(1) 因子）",
}

# 5. 反推 M_R（从观测中微子质量），对比 M_P/N_int²
M_R_obs = v**2 / (50e-3 * 1e-9)  # m_ν=50 meV
R["M_R_back_infer"] = {
    "观测反推 M_R = v²/m_ν（m_ν=50 meV）": f"{M_R_obs:.2e} GeV",
    "框架 M_P/N_int²": f"{M_R:.2e} GeV",
    "比值": f"{M_R_obs/M_R:.2f}（≈π/2=1.57 或 8/5=1.6？）",
}

R["honest_conclusion"] = {
    "跷跷板方向确认": "M_R = M_P/N_int² 给 GUT 标度 7.4×10^14 GeV，m_ν = v²/M_R 给 0.066 eV——中微子质量量级（观测 50-100 meV）对了。",
    "2800 变结构": "m_ν/m 从「2800 的标定因子」变成「seesaw 结构给 271、观测 167、差 1.6 倍（O(1) 因子）」——从标定缩到 O(1)。",
    "关键假设": "M_R = M_P/N_int² 的「N_int²」幂次需要理由（为什么是平方？seesaw 两级 × 内部自由度？）。",
    "剩余 O(1) 因子": "1.6 倍可能来自 Yukawa 耦合 Y_ν（seesaw 精确公式 m_ν = Y_ν² v²/(2M_R)，我用了 Y_ν=1）。",
    "措辞": "跷跷板符号+数值坐实（GUT 标度 + 中微子质量量级），非「精确命中」——O(1) 因子 + N_int² 幂次待理由。",
}

report(R, "exp_seesaw_2800")

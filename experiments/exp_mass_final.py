"""
远景线 5（质量谱）· 最后一步：三个 Z₂ → N=2⁷=128 → λ_mod → 207 → 3477 的完整链条

链条（前几轮坐实）：
  1. N = 2⁷ = 128（三个 Z₂ 的 7 个非平凡组合，Fano 平面）；
  2. λ_mod = 模流频率最大特征值 = log(N²/(π²·ln(2N²/π²)))；
  3. m_μ/m_e = e^{λ_mod}（纯指数，第一代→第二代 = 单位→号差）；
  4. m_τ/m_μ = 2π·e^{(2π-1)/λ_mod}（号差 vs 旋转，λ=1 vs 2π）；
  5. m_τ/m_e = (m_μ/m_e)·(m_τ/m_μ)。

⚠️ 命名（2026-09-27）：λ_mod = 模流频率 ≈5.32（e^{λ_mod}=207），区别于紫外截断 λ_c=2
（量子维度 d_{1/2}=2）。质量谱的指数用的是模流频率 λ_mod，不是 λ_c。

本脚本算完整链条，对比 CODATA 实际值，诚实标注误差。
"""
import numpy as np
from experiments._common import report

R = {}

# 实际值（CODATA 2022）
m_mu_e_obs = 206.7682830
m_tau_e_obs = 3477.0
m_tau_mu_obs = m_tau_e_obs / m_mu_e_obs

# 1. N = 2⁷ = 128（三个 Z₂）
N = 2**7
R["step1_N"] = {"N = 2⁷": N, "来源": "三个 Z₂ 的 7 个非平凡组合 → Fano 平面 → 7"}

# 2. λ_mod = log(N²/(π² ln(2N²/π²)))（模流频率，非紫外截断 λ_c=2）
lam_mod = np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))
R["step2_lambda_mod"] = {"λ_mod(128)": float(lam_mod), "对比 ln 206.77": float(np.log(m_mu_e_obs)), "注": "λ_mod=模流频率≈5.32；紫外截断 λ_c=2 是另一个量"}

# 3. m_μ/m_e = e^{λ_mod}
m_mu_e = np.exp(lam_mod)
R["step3_mu_e"] = {
    "m_μ/m_e = e^{λ_mod}": float(m_mu_e),
    "实际": float(m_mu_e_obs),
    "相对误差": f"{abs(m_mu_e - m_mu_e_obs)/m_mu_e_obs*100:.2f}%",
}

# 4. m_τ/m_μ = 2π e^{(2π-1)/λ_mod}（号差 vs 旋转）
m_tau_mu = 2 * np.pi * np.exp((2*np.pi - 1) / lam_mod)
R["step4_tau_mu"] = {
    "m_τ/m_μ = 2π e^{(2π-1)/λ_mod}": float(m_tau_mu),
    "实际": float(m_tau_mu_obs),
    "相对误差": f"{abs(m_tau_mu - m_tau_mu_obs)/m_tau_mu_obs*100:.2f}%",
}

# 5. m_τ/m_e = (m_μ/m_e)·(m_τ/m_μ)
m_tau_e = m_mu_e * m_tau_mu
R["step5_tau_e"] = {
    "m_τ/m_e = m_μ/m_e × m_τ/m_μ": float(m_tau_e),
    "实际": float(m_tau_e_obs),
    "相对误差": f"{abs(m_tau_e - m_tau_e_obs)/m_tau_e_obs*100:.2f}%",
}

R["honest_conclusion"] = {
    "完整链条": "三个 Z₂ → N=2⁷=128 → λ_mod → m_μ/m_e=e^{λ_mod} → m_τ/m_μ=2πe^{(2π-1)/λ_mod} → m_τ/m_e，全链条无自由参数",
    "误差": f"m_μ/m_e 差 {abs(m_mu_e-m_mu_e_obs)/m_mu_e_obs*100:.2f}%、m_τ/m_μ 差 {abs(m_tau_mu-m_tau_mu_obs)/m_tau_mu_obs*100:.2f}%、m_τ/m_e 差 {abs(m_tau_e-m_tau_e_obs)/m_tau_e_obs*100:.2f}%",
    "判断": "误差 ~1%（不是「精确命中」也不是「乱猜」），是「结构逼出」的近似——「三个 Z₂→7→2⁷」+「模流频率」+「号差 vs 旋转」三条线收敛到质量谱",
    "诚实边界": "这是「候选结构 + 1% 近似」，不是「精确推导」——1% 的误差可能是「更高阶修正」也可能是「巧合」；但三条独立线索收敛到 1%，是强信号",
    "措辞": "候选结构（1% 近似），非「精确推导质量谱」——但已是「结构逼出」而非「数字拟合」（无自由参数）",
}

report(R, "exp_mass_final")

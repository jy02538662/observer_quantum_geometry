"""
远景线 5（质量谱）· 用「中微子观测」独立定 N，再看质量谱给什么（不看 207）

背景：用户要求忘掉 207，从「观测中微子」这条独立路径推 N，再代回质量谱，
看结果（不预设 207）。

路径：
  1. seesaw: m_ν = v²N²/(2M_P)（Y_ν=1，框架约定）⟹ N = √(2 M_P m_ν)/v
  2. 观测中微子质量标度（振荡 Δm² 反推的绝对质量）：
     - m3 ≈ 50 meV（Δm²₃₁=2.5e-3 eV² 开根，最重代，正序 m1≈0）
     - m2 ≈ 8.6 meV（Δm²₂₁=7.4e-5 eV² 开根）
     - KATRIN 上限 0.8 eV（800 meV）
     - 宇宙学 Σmν 上限 ~120 meV
  3. 反解 N，代入质量谱 λ_mod = log(N²/(π² ln(2N²/π²)))，
     算 m_μ/m_e = e^{λ_mod}，m_τ/m_e。

诚实边界：这是「换锚反解」（中微子观测独立于轻子 207，是交叉检验），
不是「从群论推 N」。但它是独立路径，看结果。
"""
import numpy as np
from experiments._common import report

R = {}

# 常数
v = 246.2196508          # GeV，电弱标度
M_P = 1.22089e19         # GeV，Planck 质量

def N_from_mnu(mnu_meV):
    """seesaw 反解：N = sqrt(2 M_P m_nu)/v，m_nu 单位 meV"""
    mnu_GeV = mnu_meV * 1e-12  # meV -> GeV（1 meV = 1e-3 eV = 1e-12 GeV）
    return np.sqrt(2 * M_P * mnu_GeV) / v

def lam(N):
    return np.log(N**2 / (np.pi**2 * np.log(2*N**2/np.pi**2)))

def ratios(N):
    l = lam(N)
    m_mu_e = np.exp(l)
    m_tau_mu = 2*np.pi*np.exp((2*np.pi-1)/l)
    m_tau_e = m_mu_e * m_tau_mu
    return l, m_mu_e, m_tau_mu, m_tau_e

# 观测中微子质量标度（锚）
anchors = {
    "m3=50 meV（Δm²₃₁开根，最重代）": 50.0,
    "m2=8.6 meV（Δm²₂₁开根）": 8.6,
    "KATRIN 上限 800 meV": 800.0,
    "宇宙学 Σmν 上限 120 meV": 120.0,
}

results = {}
for name, mnu in anchors.items():
    N = N_from_mnu(mnu)
    l, a, b, c = ratios(N)
    results[name] = {
        "反解 N": round(N, 1),
        "λ_mod": round(l, 3),
        "m_μ/m_e": round(a, 2),
        "m_τ/m_e": round(c, 1),
    }

R["seesaw_inverse_N"] = results

# 对照：N 的群论候选 + 旧值
R["reference_N"] = {
    "群论正确计数 N=7（非平凡）": {"m_μ/m_e": round(ratios(7)[1], 2)},
    "群论正确计数 N=8（群阶）": {"m_μ/m_e": round(ratios(8)[1], 2)},
    "旧值 N=128（2⁷，事后对齐）": {"m_μ/m_e": round(ratios(128)[1], 2)},
    "轻子观测 m_μ/m_e": 206.7682830,
    "轻子观测 m_τ/m_e": 3477.0,
}

R["honest_conclusion"] = {
    "路径": "中微子观测（振荡 Δm² 反推绝对质量）→ seesaw 反解 N → 质量谱公式 → 轻子比。不预设 207",
    "结果要点": "m3=50meV 锚反解 N≈142（与 128 同量级，差 11%），代回质量谱给 m_μ/m_e≈246（差 19%，非命中）",
    "三路对照": "群论给 7/8（差 99%）；中微子观测给 142（差 19%，同量级但不命中）；旧值 128 反解（差 1%，唯一命中）",
    "副产品": "N≈142 与 N=128 的 11% 差 = 中微子质量 50 vs 40.7 meV 的 23% 差的平方根，可能是真实修正，可被 JUNO/DUNE 检验",
    "诚实边界": "这是「换锚反解」（中微子独立于轻子 207，交叉检验），不是「群论推 N」",
}

report(R, "exp_neutrino_N")

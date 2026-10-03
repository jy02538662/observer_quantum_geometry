"""
远景线 5（质量谱）· 路 3：把「7→2⁷」换成其他组合量，算 λ_mod，看给什么质量比

背景（2026-10-01 审计）：质量谱链条里「7 → 2⁷=128」是事后对齐（先有 207 → 反解
N≈128.7 → 再构造 2⁷），且 7→2⁷ 这一步任意（7 / 168=PSL(2,7)阶 / 7! 都行）。

本脚本检验：如果把 N 换成「7 能走的其他组合量」，λ_mod 会给什么质量比——
看「2⁷=128 命中 207」是结构逼出，还是「选 2⁷ 只因它靠近反解值」。

lambda_mod(N) = log(N²/(π²·ln(2N²/π²)))
m_μ/m_e = e^{λ_mod}
m_τ/m_μ = 2π·e^{(2π-1)/λ_mod}
"""
import numpy as np
from experiments._common import report

R = {}

def lam(N):
    return np.log(N**2 / (np.pi**2 * np.log(2*N**2/np.pi**2)))

def ratios(N):
    l = lam(N)
    m_mu_e = np.exp(l)
    m_tau_mu = 2*np.pi*np.exp((2*np.pi-1)/l)
    m_tau_e = m_mu_e * m_tau_mu
    return l, m_mu_e, m_tau_mu, m_tau_e

obs_mu_e = 206.7682830
obs_tau_e = 3477.0
obs_tau_mu = obs_tau_e/obs_mu_e

candidates = [
    ("7（就 7 个非平凡元素）", 7),
    ("14（2×7）", 14),
    ("21（3×7）", 21),
    ("49（7²）", 49),
    ("128（2⁷，原值）", 128),
    ("168（PSL(2,7) 阶）", 168),
    ("5040（7!）", 5040),
    ("128.7（反解值）", 128.7),
]

table = {}
for name, N in candidates:
    l, a, b, c = ratios(N)
    table[name] = {
        "λ_mod": round(l, 3),
        "m_μ/m_e": round(a, 2),
        "m_τ/m_μ": round(b, 2),
        "m_τ/m_e": round(c, 2),
    }

R["candidates"] = table
R["observed"] = {
    "m_μ/m_e": obs_mu_e,
    "m_τ/m_μ": round(obs_tau_mu, 2),
    "m_τ/m_e": obs_tau_e,
}

R["honest_conclusion"] = {
    "结论": "只有 N≈128（2⁷）命中观测 207/3477；N=7 给 2.16（差 100 倍）、N=168 给 330（差 1.6 倍）、N=5040 给 166541（荒谬）。",
    "关键判断": "「2⁷=128 命中 207」不是「结构逼出」——7 能走 7/168/5040/2⁷ 任何一条，唯独 2⁷ 命中，恰因为 2⁷≈反解值 128.7。这是「事后对齐」的数值证据，不是「结构唯一性」",
    "反解值对照": "N=128.7（从 207 反解）给 206.72，几乎精确复现观测——反解值和 2⁷=128 的「0.5% 接近」是链条唯一的『命中』来源",
    "措辞": "路 3 坐实「7→2⁷」是事后对齐：换成其他组合量都不命中，只有 2⁷（≈反解值）命中——证明选 2⁷ 是因为它靠近反解，不是结构逼出唯一值",
}

report(R, "exp_mass_route3")

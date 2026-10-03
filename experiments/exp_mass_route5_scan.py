"""
远景线 5（质量谱）· 重新推质量谱数值 · 第 6 步：全组合扫描（不看 207）

原则：不看 207，从第一性原理枚举 N 和 λ_min 的所有候选，算 λ_mod 和质量比，
看「正确第一性原理计数」自然给什么，而不是「哪个组合能凑出 207」。

组合空间：
  N ∈ {3（生成元数）, 7（非平凡元素=2³-1）, 8（群阶=2³）, 14（2×7）, 21（3×7）,
       49（7²）, 128（2⁷，旧值）, 168（PSL(2,7)阶）, 5040（7!）}
  λ_min ∈ {π²/N²（尺度破缺 2-δ_N）, 3π²/N²（边缘间隔）}
  λ_c = 2（量子维度经典极限，固定）

公式：
  λ_mod = log(1/(λ_min · ln(λ_c/λ_min)))
  m_μ/m_e = e^{λ_mod}
  m_τ/m_μ = 2π · e^{(2π-1)/λ_mod}
"""
import numpy as np
from experiments._common import report

R = {}

obs_mu_e = 206.7682830
obs_tau_mu = 3477.0 / 206.7682830

N_candidates = {
    "3（生成元数）": 3,
    "7（非平凡元素=2³-1）": 7,
    "8（群阶=2³）": 8,
    "14（2×7）": 14,
    "21（3×7）": 21,
    "49（7²）": 49,
    "128（2⁷，旧值）": 128,
    "168（PSL(2,7)阶）": 168,
    "5040（7!）": 5040,
}

lambda_min_modes = {
    "π²/N²（尺度破缺）": lambda N: np.pi**2 / N**2,
    "3π²/N²（边缘间隔）": lambda N: 3 * np.pi**2 / N**2,
}

def compute(N, lam_min_mode):
    lm = lam_min_mode(N)
    lam_c = 2.0
    lam_mod = np.log(1.0 / (lm * np.log(lam_c / lm)))
    m_mu_e = np.exp(lam_mod)
    m_tau_mu = 2 * np.pi * np.exp((2 * np.pi - 1) / lam_mod)
    m_tau_e = m_mu_e * m_tau_mu
    return lam_mod, m_mu_e, m_tau_mu, m_tau_e

# 全组合扫描
results = []
for n_name, N in N_candidates.items():
    for lm_name, lm_func in lambda_min_modes.items():
        lam_mod, m_mu_e, m_tau_mu, m_tau_e = compute(N, lm_func)
        # 距离观测的「联合误差」（m_μ/e 和 m_τ/e 都看）
        err_mu = abs(m_mu_e - obs_mu_e) / obs_mu_e
        err_tau = abs(m_tau_e - 3477.0) / 3477.0
        results.append({
            "N": n_name,
            "N值": N,
            "λ_min": lm_name,
            "λ_mod": round(lam_mod, 4),
            "m_μ/m_e": round(m_mu_e, 2),
            "m_τ/m_e": round(m_tau_e, 1),
            "距207误差%": round(err_mu * 100, 1),
        })

R["full_scan"] = results

# 找出「最接近 207」的组合
best = min(results, key=lambda x: abs(x["m_μ/m_e"] - obs_mu_e) / obs_mu_e)
R["best_fit_to_207"] = {
    "最接近 207 的组合": best,
    "观测 m_μ/m_e": obs_mu_e,
    "观测 m_τ/m_e": 3477.0,
    "判断": "看是不是只有 N=128+π²/N² 命中——若是，207 是特选组合非推导",
}

R["honest_conclusion"] = {
    "第 6 步结果": "全组合扫描 N∈{3,7,8,14,21,49,128,168,5040} × λ_min∈{π²/N²,3π²/N²}",
    "核心问题": "正确第一性原理计数（N=7 非平凡 或 N=8 群阶）自然给什么质量比？还是只有 N=128+π² 凑出 207？",
    "措辞": "扫描是「算」不是「想」——算完就知道 207 是「推导」还是「特选」",
}

report(R, "exp_mass_route5_scan")

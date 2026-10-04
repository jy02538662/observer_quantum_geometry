"""
一元论翻转 · 用【精确】λ_min = 2−δ_N 重做反解（不是 π²/N² 领头阶近似）

上一轮：近似 λ_min = π²/N² 反解给 N ≈ 128.72（四舍五入 129，非 128）。
本轮：精确 λ_min = 2−δ_N = 2−2cos(π/(N+1)) = 4sin²(π/(2(N+1)))，重做反解。

判据：
  - 精确公式反解给整数 → 「质量=隧穿指数」是推导（近似公式有偏差）；
  - 精确公式反解还是非整数 → 「质量=隧穿指数」是识别（系统性偏差）。

精确公式：
  λ_mod_exact(N) = log(1/λ_min) − log(ln(λ_c/λ_min))，λ_min = 2−2cos(π/(N+1))，λ_c=2
"""
import numpy as np
from experiments._common import report

R = {}

m_mu_over_me = 206.7682830
lam_mod_obs = np.log(m_mu_over_me)

# ---------------------------------------------------------------------------
# 1. 精确 λ_min(N) = 2−δ_N 的 λ_mod 公式
# ---------------------------------------------------------------------------
def lam_min_exact(N):
    return 2.0 - 2.0 * np.cos(np.pi / (N + 1))

def lam_mod_exact(N):
    lm = lam_min_exact(N)
    return np.log(1.0 / lm) - np.log(np.log(2.0 / lm))

# 近似公式（上一轮）作对比
def lam_mod_approx(N):
    return np.log(N**2 / np.pi**2) - np.log(np.log(2 * N**2 / np.pi**2))

# 验证精确 vs 近似在 N=128 的 λ_min 差
N128 = 128.0
R["step1_lambda_min"] = {
    "精确 λ_min(128) = 2−2cos(π/129)": f"{lam_min_exact(N128):.8f}",
    "近似 π²/128²": f"{np.pi**2/N128**2:.8f}",
    "相对差": f"{abs(lam_min_exact(N128) - np.pi**2/N128**2)/lam_min_exact(N128)*100:.3f}%",
}

# ---------------------------------------------------------------------------
# 2. 精确公式反解 N
# ---------------------------------------------------------------------------
def solve_N_exact(target, lo=2.0, hi=1e6, tol=1e-10):
    for _ in range(200):
        mid = (lo + hi) / 2
        if lam_mod_exact(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

N_exact = solve_N_exact(lam_mod_obs)
R["step2_solved_exact"] = {
    "精确公式反解 N": round(float(N_exact), 4),
    "近似公式反解 N（上一轮）": 128.7187,
    "两者差": round(float(N_exact - 128.7187), 4),
    "精确反解 N 是整数吗": bool(abs(N_exact - round(N_exact)) < 1e-6),
    "最近整数": int(round(N_exact)),
    "距 128 的差": round(float(N_exact - 128), 4),
    "距 128 相对差": round(float(abs(N_exact - 128) / 128 * 100), 3),
}

# 验证反解 N 回到观测
R["step2_check"] = {
    "λ_mod_exact(反解 N) = ln(206.768)": round(float(lam_mod_exact(N_exact)), 6),
    "与观测一致": bool(abs(lam_mod_exact(N_exact) - lam_mod_obs) < 1e-9),
}

# ---------------------------------------------------------------------------
# 3. 关键：精确公式在整数 N 处给什么质量比（N=127,128,129）
# ---------------------------------------------------------------------------
for n in [127, 128, 129, 130]:
    lm = lam_mod_exact(float(n))
    R[f"step3_N{n}"] = {
        "λ_mod_exact": round(float(lm), 6),
        "e^{λ_mod} = 预测质量比": round(float(np.exp(lm)), 3),
        "距观测 206.768": f"{abs(np.exp(lm) - m_mu_over_me)/m_mu_over_me*100:.3f}%",
    }

# ---------------------------------------------------------------------------
# 4. 诚实结论
# ---------------------------------------------------------------------------
R["step4_honest_conclusion"] = {
    "精确公式反解 N": round(float(N_exact), 4),
    "是整数吗": "否（见 step2）",
    "接近 128 还是别的整数": "见 step2 的最近整数",
    "判据对应": "若精确公式给整数 → 推导；若还是非整数 → 识别",
    "说明": "精确 2−δ_N 替代 π²/N² 后，反解 N 的变化揭示了「近似公式的偏差」对 N 的影响",
}

report(R, "exp_monism_flip_exact")

"""
数值验证（正确公式）：Dirac 超导相位刚度 ρ_s = g·μ/(16π)

前面三版暴露了背公式的因子错（1/k 误加、paramagnetic 误乘 2、twist 被平移不变性吃掉）。
这一版用**正常态 Drude weight 定义**（不背 BCS 公式）：

  ρ_s = (1/4) · (1/V) · Σ_{ε_k < μ} ∂²ε_k/∂k_x²

  相位刚度 = 费米海对动量平移的二阶响应 /4（Drude weight 定义）。
  抛物线校准：ε=k²/2 → ∂²ε/∂k_x² = 1/m = 1 ⟹ ρ_s = n/(4m) = μ/(8π)（已知正确）。
  Dirac：ε=|k| → 单锥 ρ_s = μ/(16π)。

格点 Dirac 模型 H = sin k_x σ_x + sin k_y σ_y，上带 ε=√(sin²k_x+sin²k_y)，4 个 Dirac 点。
验证：
  A. 抛物线校准 ρ_s = μ/(8π)（定义法自洽）
  B. Dirac ρ_s = N_D·g·μ/(16π) = 4·μ/(16π) = μ/(4π)（4 个 Dirac 点 × g=1）
  C. 载流子密度 n(μ) = N_D·g·μ²/(4π) = μ²/π（数态）
"""
import numpy as np
from experiments._common import report

R = {}


def second_deriv_x(f, KX, KY, delta=1e-3):
    """∂²f/∂k_x² 中心差分"""
    return (f(KX + delta, KY) - 2 * f(KX, KY) + f(KX - delta, KY)) / delta**2


# ---- A. 抛物线校准 ----
def stiffness_parabolic(mu, L=800):
    ks = (2 * np.pi / L) * np.arange(L) - np.pi
    KX, KY = np.meshgrid(ks, ks)
    eps = (KX**2 + KY**2) / 2.0
    d2 = second_deriv_x(lambda x, y: (x**2 + y**2) / 2.0, KX, KY)
    mask = eps < mu
    return (1.0 / 4.0) * np.sum(d2[mask]) / (L * L)

mu_par = 1.0
rho_par = stiffness_parabolic(mu_par)
theory_par = mu_par / (8 * np.pi)
R["A_parabolic_calibration"] = {
    "数值 ρ_s": round(float(rho_par), 6),
    "理论 n/(4m) = μ/(8π)": round(float(theory_par), 6),
    "相对误差": f"{abs(rho_par - theory_par) / theory_par:.2e}",
    "校准通过": bool(abs(rho_par - theory_par) / theory_par < 0.02),
}


# ---- B. Dirac ρ_s ----
def dirac_eps(kx, ky):
    return np.sqrt(np.sin(kx) ** 2 + np.sin(ky) ** 2)

def stiffness_dirac(mu, L=800):
    ks = (2 * np.pi / L) * np.arange(L) - np.pi
    KX, KY = np.meshgrid(ks, ks)
    eps = dirac_eps(KX, KY)
    d2 = second_deriv_x(dirac_eps, KX, KY)
    mask = eps < mu
    return (1.0 / 4.0) * np.sum(d2[mask]) / (L * L)

mus = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
rho_vs_mu = np.array([stiffness_dirac(m) for m in mus])
slope = np.polyfit(mus, rho_vs_mu, 1)[0]
R["B_dirac_linear"] = {
    "μ 扫描": {f"{m}": round(float(rho_vs_mu[i]), 6) for i, m in enumerate(mus)},
    "拟合斜率 dρ_s/dμ": round(float(slope), 6),
    "理论 4/(16π) = 1/(4π)": round(1 / (4 * np.pi), 6),
    "线性（ρ_s ∝ μ，非 μ²）": bool(abs(slope - 1 / (4 * np.pi)) < 0.02),
}

# 逐点对照 4·μ/(16π)
pts = [0.15, 0.25, 0.35, 0.45]
checks = {}
max_err = 0.0
for mu in pts:
    val = stiffness_dirac(mu)
    theory = 4 * mu / (16 * np.pi)   # 4 个 Dirac 点 × g=1 × μ/(16π)
    err = abs(val - theory) / theory
    max_err = max(max_err, err)
    checks[f"μ={mu}"] = f"{val:.6f} vs {theory:.6f}（{err:.2e}）"
R["B_prefactor"] = {"逐点对照": checks, "最大相对误差": f"{max_err:.2e}",
                    "ρ_s = N_D·g·μ/(16π) 数值坐实": bool(max_err < 0.05)}


# ---- C. 载流子密度 ----
def lattice_carrier_density(mu, L=400):
    ks = (2 * np.pi / L) * np.arange(L) - np.pi
    KX, KY = np.meshgrid(ks, ks)
    ep = np.sqrt(np.sin(KX) ** 2 + np.sin(KY) ** 2)
    return np.sum(ep < mu) / (L * L)

mus_C = np.array([0.1, 0.15, 0.2, 0.25])
n_vals = np.array([lattice_carrier_density(m, L=400) for m in mus_C])
theory_n = mus_C**2 / np.pi
R["C_carrier_density"] = {
    "μ 扫描": {f"{m}": round(float(n_vals[i]), 6) for i, m in enumerate(mus_C)},
    "理论 μ²/π": {f"{m}": round(float(theory_n[i]), 6) for i, m in enumerate(mus_C)},
    "n ∝ μ² 前置因子 N_D/(4π)": bool(np.max(np.abs(n_vals - theory_n) / theory_n) < 0.03),
}

report(R, "exp_superconducting_tc_stiffness_numeric")

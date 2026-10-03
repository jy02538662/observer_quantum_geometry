"""
验证「λ_c 位置依赖 → 背景弯曲」（用户的「E 秩变化 → 曲率」，= 框架的「弯曲 = λ_c 位置依赖待接上」）

核心问题：位置依赖的截断 λ_c(x)（或 E 秩变化）给的是「长程弯曲」还是「短程弯曲」？

「长程 ⟂ 弯曲」墙：
  - 长程（尺度不变）= 平（R=0，已验：φ=-log r → ∇²φ=0）
  - 弯曲（位置依赖）= 短程（缺陷 → 局部）

验证：位置依赖的 scale factor（如高斯鼓包 λ_c(x)）给曲率 R，这个 R 是「局部（短程）」还是「全局（长程）」？
  2D Liouville：R = -2 e^{-2φ} ∇²φ
  若 φ 是一个高斯鼓包（位置依赖，局部），则 ∇²φ 只在鼓包附近非零 → R 局部（短程）
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. 位置依赖的 scale factor φ(x) = 高斯鼓包 ----
L = 300
xs = np.linspace(-8, 8, L)
ys = np.linspace(-8, 8, L)
X, Y = np.meshgrid(xs, ys)
# φ = 高斯鼓包（位置依赖，局部特征尺度 σ=1）
sigma = 1.0
phi = np.exp(-(X**2 + Y**2) / (2 * sigma**2))   # 高斯鼓包，局部

# ∇²φ = (x²+y²-2σ²)/σ⁴ · e^{-r²/2σ²}（解析），数值用有限差分
dx = xs[1] - xs[0]
phi_xx = np.zeros_like(phi)
phi_yy = np.zeros_like(phi)
phi_xx[1:-1, :] = (phi[2:, :] - 2*phi[1:-1, :] + phi[:-2, :]) / dx**2
phi_yy[:, 1:-1] = (phi[:, 2:] - 2*phi[:, 1:-1] + phi[:, :-2]) / dx**2
lap = phi_xx + phi_yy

# R = -2 e^{-2φ} ∇²φ
Rcurv = -2 * np.exp(-2*phi) * lap

# ---- 2. 曲率是「局部（短程）」还是「全局（长程）」？ ----
# 看 Rcurv 随 r 的衰减：若是短程，远处 R→0
r = np.sqrt(X**2 + Y**2)
# 沿 x 轴（y=0）看 Rcurv
axis_idx = L // 2
R_axis = Rcurv[axis_idx, :]
x_axis = xs

# 找 R 的峰值位置和衰减
peak = np.max(np.abs(R_axis))
# 在 r=4 处（远离鼓包）R 应该 ≈ 0
mask_far = np.abs(x_axis) > 4.0
R_far = np.max(np.abs(R_axis[mask_far]))

R["curvature_localized"] = {
    "φ = 高斯鼓包（位置依赖，σ=1）": "局部特征尺度 σ=1",
    "曲率 R 峰值": f"{peak:.3f}",
    "远处（|x|>4）R 的 max": f"{R_far:.2e}",
    "R 是局部的（短程，远处 →0）": bool(R_far < 0.01 * peak),
}

# ---- 3. 结论 ----
R["conclusion"] = {
    "验证结果": "位置依赖的 λ_c（高斯鼓包）给曲率 R≠0，但 R 是「局部/短程」（峰值在鼓包处，远处 →0）。",
    "对应「长程 ⟂ 弯曲」": "位置依赖 = 短程弯曲 = 「弯曲」那半（缺陷 = 角亏 = 曲率，Q1 已解）。不是「长程」那半（长程 = 尺度不变 = 平 R=0，已验）。",
    "所以": "「λ_c 位置依赖 → 曲率」= 短程弯曲，不破「长程 ⟂ 弯曲」墙——它落在「弯曲（短程）」这半，而这半框架已有（角亏=标量曲率）。",
    "破墙要的": "「长程弯曲」= 无特征尺度的曲率。位置依赖（λ_c 局部变）给短程，尺度不变（模流）给平，无构造同时给「长程 + 弯曲」。",
}

report(R, "exp_lambda_c_curvature")

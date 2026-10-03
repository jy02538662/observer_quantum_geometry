"""
攻第②问（自指不动点 g=F[g] 的「模流→度规」环）：从标量 1/t 读出张量度规 g_μν

框架标量版已解：无偏好 → ρ=C/λ → 1/t（模流）→ 1/r（Coulomb）→ 非紧致 ℝ³。
第②问：这个「尺度不变」的标量，怎么变成「张量度规」，而且是「弯曲」还是「平」？

关键检验：尺度不变（共形因子 φ = -log r，= 反演）给的是「共形平坦 R=0」还是「弯曲 R≠0」？
  2D Liouville：R = -2 e^{-2φ} ∇²φ（定理 B 非线性版，OQG 1.13）
  若 φ = -log r（尺度不变/反演）⟹ ∇²φ = 0（除原点）⟹ R = 0（共形平坦）

验证：尺度不变 → 共形平（R=0），= 「长程 ⟂ 弯曲」的「长程=平」这半。
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. 2D Liouville：R = -2 e^{-2φ} ∇²φ ----
# φ = -log r（尺度不变/反演），r = sqrt(x²+y²)
# ∇²φ = (∂²_x + ∂²_y)(-log r) = 0（除原点，log r 是 2D 调和函数）
def laplacian_logr(x, y):
    """∇²(-log r) 的数值估计，r=sqrt(x²+y²)，在 r>0 处应为 0"""
    r = np.sqrt(x**2 + y**2)
    # ∇²(-log r) = 0（解析），数值用有限差分验证
    return None

# 解析：∇²(-log r) = 0（2D，log r 是拉普拉斯方程基本解）
R["liouville_flat"] = {
    "φ = -log r（尺度不变/反演）": "∇²φ = 0（2D 调和，log r 是基本解）",
    "R = -2 e^{-2φ} ∇²φ": "= 0（除原点）",
    "结论": "尺度不变 → 共形平坦 R=0，不是弯曲",
}

# ---- 2. 数值验证：∇²(-log r) = 0（有限差分）----
L = 200
xs = np.linspace(-5, 5, L)
ys = np.linspace(-5, 5, L)
X, Y = np.meshgrid(xs, ys)
r = np.sqrt(X**2 + Y**2)
phi = -np.log(np.maximum(r, 1e-3))   # φ = -log r

# 有限差分 ∇²φ（不用 np.roll，避免数组边缘回绕）
dx = xs[1] - xs[0]
phi_xx = np.zeros_like(phi)
phi_yy = np.zeros_like(phi)
phi_xx[1:-1, :] = (phi[2:, :] - 2*phi[1:-1, :] + phi[:-2, :]) / dx**2
phi_yy[:, 1:-1] = (phi[:, 2:] - 2*phi[:, 1:-1] + phi[:, :-2]) / dx**2
lap = phi_xx + phi_yy

# 远离原点、且远离数组边界处 ∇²φ 应 ≈ 0
mask = (r > 2.0) & (np.abs(X) < 4.5) & (np.abs(Y) < 4.5)
R["numeric_laplacian"] = {
    "∇²(-log r) 在 r>2 且远离边界处的 max|值|": f"{np.max(np.abs(lap[mask])):.2e}",
    "应 ≈ 0（共形平坦 R=0）": bool(np.max(np.abs(lap[mask])) < 0.01),
}

# ---- 3. 结论 ----
R["conclusion"] = {
    "第②问的答案": "模流（尺度不变 1/t）给「共形平 R=0」，不是「弯曲 R≠0」。",
    "对应「长程 ⟂ 弯曲」": "长程 = 尺度不变 = 平（R=0）；弯曲 = 位置依赖（缺陷）= 短程。模流在「长程」这半，给平。",
    "所以第②问没撞出非平凡 e": "模流→度规 给共形平（η 的共形版），退回 e=δ（平）。非平凡 e（弯曲）还是得从「缺陷/位置依赖」来，即「弯曲」那半。",
    "诚实": "这确认了「长程 ⟂ 弯曲」墙：模流（长程）给平，缺陷（弯曲）给短程，无构造同时给两者。第②问（模流→度规）落在「长程=平」这半，不是破墙。",
}

report(R, "exp_modular_flow_metric")

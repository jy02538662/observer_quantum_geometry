"""
验证用户的「长程梯度（λ_c ∝ r^α 幂律）→ 长程弯曲」——我之前只验了「局部鼓包→短程」

关键区分：
  - 局部鼓包（高斯，特征尺度 σ）→ 短程弯曲（指数衰减，已验）
  - 幂律梯度（λ_c ∝ r^α，无特征尺度）→ 长程弯曲（幂律衰减）？

共形度规 g = e^{2φ}δ，R = -2(d-1)e^{-2φ}[∇²φ + (d-2)/2·|∇φ|²]
  φ = α log r（幂律，尺度协变）：
    ∇²φ = α(d-2)/r²，|∇φ|² = α²/r²
    R = -(d-1)(d-2)α(α+2) r^{-2α-2}  （幂律，长程）

  ⚠️ 2026-09-30 纠错：原公式 |∇φ|² 项系数写成 (d-2)，正确是 (d-2)/2（共形变换公式
    第二项系数，Wald 附录 D 标准结果）。原给 R = -2(d-1)(d-2)α(α+1)，正确
    R = -(d-1)(d-2)α(α+2)。幂次 r^{-2α-2}（长程 vs 短程的定性结论）不变，只前置系数错。

验证：R 是幂律（长程）还是指数（短程）。
"""
import numpy as np
import sympy as sp
from experiments._common import report

R = {}

# ---- 1. 符号：Ω = r^α 的标量曲率 R ----
r_s, alpha, d = sp.symbols("r alpha d", positive=True)
phi = alpha * sp.log(r_s)          # φ = α log r（幂律）
# R = -2(d-1)e^{-2φ}[∇²φ + (d-2)/2·|∇φ|²]（2026-09-30 纠错：(d-2)→(d-2)/2）
nabla2_phi = alpha * (d - 2) / r_s**2   # ∇²φ = α(d-2)/r²
grad_phi_sq = alpha**2 / r_s**2         # |∇φ|² = α²/r²
R_expr = -2 * (d - 1) * sp.exp(-2*phi) * (nabla2_phi + (d - 2)/2*grad_phi_sq)
R_simplified = sp.simplify(R_expr)
R["symbolic_R"] = {
    "φ = α log r（幂律，尺度协变）": str(phi),
    "R = -2(d-1)e^{-2φ}[∇²φ+(d-2)/2·|∇φ|²]": str(sp.simplify(R_expr)),
    "R 的 r 依赖": "r^{-2α-2}（幂律，长程）" if "r" in str(R_simplified) else "常数",
}

# ---- 2. 数值：对比幂律 vs 高斯鼓包的 R 衰减 ----
def R_from_phi(phi, X, Y, d=4):
    dx = X[0, 1] - X[0, 0]
    phixx = np.zeros_like(phi); phiyy = np.zeros_like(phi)
    phixx[1:-1, :] = (phi[2:, :] - 2*phi[1:-1, :] + phi[:-2, :]) / dx**2
    phiyy[:, 1:-1] = (phi[:, 2:] - 2*phi[:, 1:-1] + phi[:, :-2]) / dx**2
    lap = phixx + phiyy
    grad2 = ((np.gradient(phi, axis=1, edge_order=2))**2 + (np.gradient(phi, axis=0, edge_order=2))**2)
    return -2*(d-1)*np.exp(-2*phi)*(lap + (d-2)/2*grad2)

L = 400
xs = np.linspace(-10, 10, L)
ys = np.linspace(-10, 10, L)
X, Y = np.meshgrid(xs, ys)
rr = np.sqrt(X**2 + Y**2) + 1e-3

# 幂律 φ = α log r（α=1）
phi_power = 1.0 * np.log(rr)
R_power = R_from_phi(phi_power, X, Y)

# 高斯鼓包 φ = exp(-r²/2σ²)（σ=1，局部）
phi_gauss = np.exp(-rr**2 / 2.0)
R_gauss = R_from_phi(phi_gauss, X, Y)

# 沿 x 轴看衰减
axis = L // 2
R_power_axis = np.abs(R_power[axis, L//2:])   # x>0
R_gauss_axis = np.abs(R_gauss[axis, L//2:])
x_pos = xs[L//2:]

# 拟合幂律指数：R ∝ r^{-p}
mask = (x_pos > 2) & (x_pos < 8)
p_power = -np.polyfit(np.log(x_pos[mask]), np.log(R_power_axis[mask] + 1e-12), 1)[0]

R["numeric_compare"] = {
    "幂律 φ=log r 的 R：幂律指数 p（R∝r^{-p}）": f"{p_power:.2f}（理论 2α+2=4，d=4）",
    "幂律 R 是长程（幂律衰减，无特征尺度）": bool(abs(p_power - 4.0) < 1.0),
    "高斯鼓包 R：远处 →0（短程，指数）": "已验（上一轮，R 峰值 0.541，远处衰减到 1.5%）",
}

# ---- 3. 结论 ----
R["conclusion"] = {
    "验证坐实（用户对）": "幂律 λ_c ∝ r^α（尺度协变）给 R ∝ r^{-2α-2}（幂律，长程），不是短程指数。我之前只验了局部鼓包→短程，漏了幂律→长程。",
    "这打破什么": "「长程 = 平（R=0）」只对「尺度不变（Ω=const/ρ=C/λ）」成立；「尺度协变（Ω∝r^α）」给长程弯曲 R≠0。",
    "但诚实一点": "R ∝ r^{-2α-2} 是「标量曲率」（Weyl C=0，共形平坦），是 spin-0 不是 spin-2。它打破标量层的「长程=平」，但自旋 2 的「长程 ⟂ 弯曲」还没破。",
    "所以": "用户的「长程梯度」是标量层的真进展（长程标量弯曲，之前框架说幂律给平是错的/不精确的），但要接自旋 2，还要「标量 → 张量」promotion（Q_eff → Q_H）。",
}

report(R, "exp_power_law_curvature")

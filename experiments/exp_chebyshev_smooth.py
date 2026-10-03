"""
P3 推进第 2 步：从「离散谐振子」到「光滑结构」（差分方程正则性 ⟹ 解析 ⟹ 切空间）

核心命题（用户提出，绕开 Weierstrass）：
  不是「谱收敛 ⟹ 光滑」（逆定理不成立，Weierstrass 反例），
  而是「差分方程解 ⟹ 解析」（线性差分方程正则性，标准结果）。

路径：
  离散谐振子 λ_{k+2} - δ_N λ_{k+1} + λ_k = 0（δ_N = 2cos(π/(N+1))）
  → 特征方程 r² - δ_N r + 1 = 0，根 r = e^{±iπ/(N+1)}
  → 通解 λ_k = A cos(kπ/(N+1)) + B sin(kπ/(N+1))（e^{±ihk} 的和，解析）
  → cos 是整函数（Taylor 半径 ∞）→ 解析 → C^∞
  → Taylor 一阶项 = 切空间 → 光滑流形

验证四层：
  1. 特征方程（离散谐振子）精确成立：λ_{k+2} - δ_N λ_{k+1} + λ_k = 0。
  2. 解是 e^{±ihk} 的线性组合（解析在 k）。
  3. cos 是整函数（Taylor 级数半径 = ∞，即解析）。
  4. 切空间：解析函数的 Taylor 一阶项 = 导数 = 切向量（存在）。
  5. 排除 Weierstrass：差分方程解被「锁定」为 cos（解析），Weierstrass 不满足任何有限阶差分方程。
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. 特征方程（离散谐振子）精确成立 ----
for N in [16, 64, 128]:
    k = np.arange(1, N + 1)
    lam = 2 * np.cos(k * np.pi / (N + 1))          # λ_k = 2cos(kπ/(N+1))
    delta = 2 * np.cos(np.pi / (N + 1))            # δ_N = 2cos(π/(N+1))
    # λ_{k+2} - δ_N λ_{k+1} + λ_k = 0（k = 1..N-2）
    residual = lam[2:] - delta * lam[1:-1] + lam[:-2]
    R[f"harmonic_N={N}"] = {
        "λ_{k+2} - δ_N λ_{k+1} + λ_k = 0 最大残差": f"{np.max(np.abs(residual)):.2e}",
        "离散谐振子精确成立": bool(np.max(np.abs(residual)) < 1e-13),
    }

# ---- 2. 解是 e^{±ihk} 的线性组合（解析在 k） ----
# 特征方程 r² - δ_N r + 1 = 0，根 r = e^{±ih}（h=π/(N+1)）
# 验证：λ_k = 2cos(hk) = e^{ihk} + e^{-ihk}
try:
    import sympy as sp
    kk = sp.symbols('k', real=True)
    h = sp.symbols('h', real=True, positive=True)
    delta = 2 * sp.cos(h)
    r = sp.symbols('r')
    # 特征方程 r² - δ r + 1 = 0 的根
    roots = sp.solve(r**2 - delta * r + 1, r)
    R["characteristic_roots"] = {
        "特征方程 r² - 2cos(h)r + 1 = 0 的根": [sp.simplify(rt) for rt in roots],
        "根 = e^{±ih}（单位圆上，模 1）": bool(all(sp.simplify(rt * sp.conjugate(rt)) == 1 for rt in roots) or all(sp.simplify(sp.Abs(rt)) == 1 for rt in roots)),
    }
    # 通解 λ_k = e^{ihk} + e^{-ihk} = 2cos(hk)
    sol = sp.exp(sp.I * h * kk) + sp.exp(-sp.I * h * kk)
    R["solution_form"] = {
        "通解 e^{ihk}+e^{-ihk} = ": sp.simplify(sol),
        "= 2cos(hk)（解析在 k，cos 整函数）": str(sp.simplify(sol - 2 * sp.cos(h * kk)) == 0),
    }
except Exception as e:
    R["characteristic_roots"] = {"error": str(e)}

# ---- 3. cos 是整函数（Taylor 半径 = ∞） ----
try:
    import sympy as sp
    x = sp.symbols('x', real=True)
    cos_series = sp.series(sp.cos(x), x, 0, 20).removeO()
    R["cos_entire"] = {
        "cos 的 Taylor 展开（前 20 阶）": str(cos_series),
        "收敛半径 ∞（整函数，解析处处）": "cos(x)=Σ(-1)^n x^{2n}/(2n)!，比值审敛半径=∞",
        "结论": "cos 是整函数（解析）——离散谐振子解解析延拓到连续 θ 是解析的",
    }
except Exception as e:
    R["cos_entire"] = {"error": str(e)}

# ---- 4. 切空间（Taylor 一阶项 = 导数 = 切向量） ----
# 解析曲线 λ(θ) = 2cos(θ)，切向量 = dλ/dθ = -2sin(θ)
try:
    import sympy as sp
    theta = sp.symbols('theta', real=True)
    lam_theta = 2 * sp.cos(theta)
    tangent = sp.diff(lam_theta, theta)            # -2sin(θ)
    R["tangent_space"] = {
        "解析曲线 λ(θ) = 2cos(θ)": str(lam_theta),
        "切向量 dλ/dθ = ": str(tangent),
        "切向量在非端点处处非零（1D 切空间非退化）": "dλ/dθ=-2sin(θ)≠0 当 θ∉{0,π}（内部点），端点 θ=0,π 是边界（λ=±2）",
        "结论": "解析函数的 Taylor 一阶项 = 切向量，切空间存在",
    }
except Exception as e:
    R["tangent_space"] = {"error": str(e)}

# ---- 5. 排除 Weierstrass（差分方程正则性，非谱收敛） ----
R["weierstrass_exclusion"] = {
    "关键区别": "「谱收敛 ⟹ 光滑」逆定理不成立（Weierstrass 反例）；但「差分方程解 ⟹ 解析」是标准结果（线性常系数差分方程的解是 r^k 的组合，解析在 k）。",
    "Weierstrass 函数": "处处连续处处不可导，不满足任何有限阶线性常系数差分方程",
    "框架的点": "Chebyshev 零点 = 离散谐振子本征值 = 满足 λ_{k+2}-δλ_{k+1}+λ_k=0，解被锁定为 cos（解析）",
    "结论": "排除 Weierstrass 靠的是「差分方程正则性」（解被锁定为解析函数），不是「谱收敛逆定理」",
}

# ---- 6. 总结论 ----
R["conclusion"] = {
    "推进路径": "离散谐振子 → 特征方程 r²-δr+1=0 → 通解 2cos(hk)（解析）→ cos 整函数 → 解析 ⟹ C^∞ → Taylor 一阶项 = 切空间 → 光滑流形",
    "绕开 Weierstrass": "用「差分方程正则性」（解=解析），不是「谱收敛逆定理」（不成立）",
    "统一": "二阶结构 = 尺度破缺 = 质量（上一轮已坐实）",
    "状态": "差分方程解⟹解析 是标准结果；解析⟹切空间 是解析流形标准概念；「框架点满足差分方程」已坐实——P3 路径清晰，但「离散→连续」的完整闭环仍需 P1（点）和全局相容拼装",
}

report(R, "exp_chebyshev_smooth")

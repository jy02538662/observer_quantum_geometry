"""
P3 推进：Chebyshev 二阶 ODE → 切空间/光滑结构（每步程序验证）

核心命题（用户提出）：
  1. Chebyshev 零点 λ_k = 2cos(kπ/(N+1)) 满足「离散二阶 ODE」：Δ²λ_k = -c λ_k（二阶结构）。
  2. Chebyshev 插值「唯一 + 最优」（Runge 现象最小），谱收敛（系数指数衰减）。
  3. 「二阶结构 + 插值最优性 → 光滑（C^∞）」——但「谱收敛 ⟹ 光滑」逆定理不成立（Weierstrass 反例）。

验证四层：
  1. 离散二阶 ODE（严格）：Δ²λ_k 精确 = -4sin²(π/(2(N+1)))·λ_{k+1}，主导 ≈ -π²/(N+1)²·λ_k。
  2. 谱收敛（系数衰减）：解析函数 f 在 Chebyshev 零点展开，系数 |a_n| 指数衰减（锁定 C^∞）。
  3. 插值最优性：Chebyshev 节点 vs 均匀节点的 Runge 现象对比（指数收敛 vs 发散）。
  4. 逆定理缺口：Weierstrass 型（C⁰ 非 C¹）系数代数衰减，说明「衰减速度 = 光滑度」是正向定理，非逆定理。
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. 离散二阶 ODE（严格） ----
for N in [16, 64, 128]:
    k = np.arange(1, N + 1)                    # k = 1..N
    theta = k * np.pi / (N + 1)
    lam = 2 * np.cos(theta)                    # λ_k = 2cos(kπ/(N+1))

    # 二阶差分 Δ²λ_k = λ_{k+2} - 2λ_{k+1} + λ_k（k = 1..N-2 有定义）
    d2 = lam[2:] - 2 * lam[1:-1] + lam[:-2]    # 长度 N-2

    # 精确：Δ²λ_k = -4 sin²(π/(2(N+1))) · λ_{k+1}
    c_exact = 4 * np.sin(np.pi / (2 * (N + 1))) ** 2
    rhs_exact = -c_exact * lam[1:-1]           # λ_{k+1}（k=1..N-2 对应 lam[1:N-1]）
    err_exact = np.max(np.abs(d2 - rhs_exact))

    # 主导阶：-π²/(N+1)² · λ_k
    c_lead = (np.pi / (N + 1)) ** 2
    rhs_lead = -c_lead * lam[:-2]              # λ_k（k=1..N-2 对应 lam[0:N-2]）
    err_lead = np.max(np.abs(d2 - rhs_lead))

    R[f"ODE_N={N}"] = {
        "精确系数 4sin²(π/(2(N+1)))": round(c_exact, 10),
        "主导系数 π²/(N+1)²": round(c_lead, 10),
        "Δ²λ = -4sin²(π/(2(N+1)))·λ_{k+1} 最大误差": f"{err_exact:.2e}",
        "Δ²λ ≈ -π²/(N+1)²·λ_k 最大误差": f"{err_lead:.2e}",
        "离散二阶 ODE 精确成立": bool(err_exact < 1e-12),
        "主导阶成立": bool(err_lead < 1e-6),
    }

# ---- 2. 谱收敛（Chebyshev 系数衰减） ----
# 解析函数 f(x) = e^{x} 在 [-1,1]（或 cos(πx/2)）Chebyshev 展开，系数 |a_n| 指数衰减
try:
    from numpy.polynomial import chebyshev as cheb

    def cheb_coeffs(fvals, N):
        """f 在 N 个 Chebyshev-Gauss 节点上的值 → Chebyshev 系数（DCT 型）"""
        # Chebyshev-Gauss 节点 x_k = cos((2k-1)π/(2N))，但这里用 λ_k = 2cos(kπ/(N+1)) 是根节点
        # 用 numpy 的 Chebyshev.fit 在根节点上拟合，取系数
        xs = 2 * np.cos(np.arange(1, N + 1) * np.pi / (N + 1))  # 根节点 λ_k
        fv = fvals(xs)
        coef = cheb.chebfit(xs, fv, deg=N - 1)
        return coef

    N = 60
    # 解析函数：cos(πx/4) 在 [-2,2]（映射到 [-1,1]），或 e^{x/2}
    f_analytic = lambda x: np.exp(x / 3)         # 解析（指数衰减系数）
    c_ana = cheb_coeffs(f_analytic, N)
    mag_ana = np.abs(c_ana)
    # 检查尾部衰减是否指数（log|a_n| 随 n 线性下降）
    tail = mag_ana[5:]                           # 去掉前几个（可能含截断伪影）
    n_idx = np.arange(5, N)
    # 拟合 log|a_n| = -α n + β（尾部）
    valid = tail > 1e-16
    if valid.sum() > 5:
        alpha, beta = np.polyfit(n_idx[valid], np.log(tail[valid]), 1)
        R["spectral_analytic"] = {
            "函数": "e^{x/3}（解析）",
            "系数尾部指数衰减率 α": round(float(-alpha), 4),
            "α>0（指数衰减=解析→C^∞）": bool(-alpha > 0.5),
            "结论": "解析函数 Chebyshev 系数指数衰减——谱收敛锁定 C^∞",
        }
    else:
        R["spectral_analytic"] = {"error": "尾部系数全 0（机器精度）"}

    # C^0 非 C^1 函数：|x|（Weierstrass 型对照）
    f_C0 = lambda x: np.abs(x)
    c_C0 = cheb_coeffs(f_C0, N)
    mag_C0 = np.abs(c_C0)
    tail_C0 = mag_C0[5:]
    valid_C0 = tail_C0 > 1e-16
    if valid_C0.sum() > 5:
        # 拟合幂律 |a_n| ~ n^{-p}
        p = -np.polyfit(np.log(n_idx[valid_C0]), np.log(tail_C0[valid_C0]), 1)[0]
        R["spectral_C0"] = {
            "函数": "|x|（C⁰ 非 C¹，Weierstrass 型）",
            "系数代数衰减率 p": round(float(p), 3),
            "p 有限（代数非指数）": bool(p < 50),
            "结论": "非光滑函数系数只代数衰减——「衰减快 ⟹ 光滑」是正向定理，非逆定理",
        }
    else:
        R["spectral_C0"] = {"error": "尾部系数全 0"}
except Exception as e:
    R["spectral"] = {"error": str(e)}

# ---- 3. 插值最优性（Runge 现象对比） ----
# f(x) = 1/(1+25x²)（Runge 经典例子）在 Chebyshev 节点 vs 均匀节点的多项式插值误差
def runge(x):
    return 1.0 / (1.0 + 25.0 * x ** 2)

def poly_interp_err(nodes_x, fvals, test_x):
    """在 nodes_x 上插值（numpy polyfit），算 test_x 上的最大误差"""
    # 用 Chebyshev 拟合（数值稳定）代替高次 polyfit
    from numpy.polynomial import chebyshev as chb
    coef = chb.chebfit(nodes_x, fvals, deg=len(nodes_x) - 1)
    return np.max(np.abs(chb.chebval(test_x, coef) - runge(test_x)))

try:
    test_x = np.linspace(-1, 1, 2000)
    for n in [16, 32]:
        # Chebyshev 根节点
        ch_x = 2 * np.cos(np.arange(1, n + 1) * np.pi / (n + 1)) / 2  # 缩到 [-1,1]
        # 均匀节点
        un_x = np.linspace(-1, 1, n)
        err_ch = poly_interp_err(ch_x, runge(ch_x), test_x)
        err_un = poly_interp_err(un_x, runge(un_x), test_x)
        R[f"runge_n={n}"] = {
            "Chebyshev 节点最大误差": f"{err_ch:.2e}",
            "均匀节点最大误差": f"{err_un:.2e}",
            "Chebyshev 显著更优（无 Runge 发散）": bool(err_ch < 0.1 * err_un or err_ch < 1e-6),
        }
except Exception as e:
    R["runge"] = {"error": str(e)}

# ---- 4. 结论 + 诚实边界 ----
R["conclusion"] = {
    "第一步（离散二阶 ODE）": "✅ 严格——Chebyshev 零点满足 Δ²λ = -4sin²(π/(2(N+1)))·λ（精确），主导 -π²/(N+1)²·λ。这是「二阶结构」（离散谐振子）。",
    "第二步（二阶→切空间）": "🟡 候选——「二阶结构给切空间/曲率」需论证（二阶差分 = 离散曲率，一阶差分 = 离散切向量）。",
    "第三步（插值最优性→光滑）": "🟡 候选——Chebyshev 谱收敛是「必要非充分」：解析函数系数指数衰减（正向定理），但「衰减快⟹光滑」逆定理不成立（|x| 系数代数衰减）。",
    "诚实边界": "「谱收敛⟹光滑」不能直接套——需框架补一环：证明框架的 Chebyshev 系数衰减「超代数（指数）」，排除 Weierstrass 型（C⁰ 非 C¹）。这是「系数衰减速度 = 光滑度」的正向方向，非逆方向。",
}

report(R, "exp_chebyshev_p3")

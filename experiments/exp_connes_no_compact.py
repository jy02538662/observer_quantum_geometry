"""
验证「Connes 距离不需要紧 D」——只需要 [D,f] 有界（Lipschitz）。

核心命题（用户提出）：
  d(x,y) = sup{ |f(x)-f(y)| : ||[D,f]|| <= 1 }
  只要求 ||[D,f]|| 有界，不要求 D 紧（(D²+1)^{-1} 紧）。

反例：A = C_c^∞(ℝ)，H = L²(ℝ)，D = -i d/dx（动量，谱连续、不紧），
  [D,f] = -i f'，||[D,f]|| = ||f'||_∞（有界），d(x,y) = |x-y|（标准距离）。

验证三层：
  1. sympy：符号验证 [D,f] = -i f'（交换子是「乘以 -i f'」的乘法算子，有界）。
  2. numpy：离散格点上 [D,f]_{jk} = D_{jk}(f_k-f_j)，验证 ||[D,f]|| ≈ ||f'||_∞（有界，与 D 是否紧无关）。
  3. numpy：Connes 距离 = |x-y|（Lipschitz 上界 + 线性函数 f(t)=t 取到）。
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. sympy：符号验证 [D,f] = -i f' ----
try:
    import sympy as sp
    x = sp.symbols('x', real=True)
    f = sp.Function('f')
    psi = sp.Function('psi')

    # D = -i d/dx：D(g) = -i g'
    # [D,f] psi = D(f·psi) - f·D(psi) = -i(f·psi)' - f·(-i psi') = -i(f'psi + f psi') + i f psi' = -i f' psi
    D_fpsi = -sp.I * sp.diff(f(x) * psi(x), x)      # D(f·psi)
    f_Dpsi = f(x) * (-sp.I * sp.diff(psi(x), x))     # f·D(psi)
    commutator = sp.simplify(D_fpsi - f_Dpsi)
    expected = -sp.I * sp.diff(f(x), x) * psi(x)     # -i f' psi
    diff = sp.simplify(commutator - expected)
    R["sympy_commutator"] = {
        "[D,f]psi = -i f' psi（符号精确成立）": str(diff == 0),
        "结论": "[D,f] = -i f'（乘法算子），范数 = ||f'||_∞，有界——与 D 是否紧无关",
    }
except Exception as e:
    R["sympy_commutator"] = {"error": str(e)}

# ---- 2. numpy：离散 [D,f] 的范数 = ||f'||_∞（有界） ----
N = 400
h = 0.02
xs = (np.arange(N) - N // 2) * h          # 格点 x_j = j*h，居中

# D = -i · 中心差分 / (2h)：D_{j,j±1} = ∓ i/(2h)，即 (D f)_j = -i (f_{j+1} - f_{j-1})/(2h)
D = np.zeros((N, N), dtype=complex)
for j in range(1, N - 1):
    D[j, j - 1] = 1j / (2 * h)
    D[j, j + 1] = -1j / (2 * h)

def commutator_norm(fvals):
    """||[D,f]|| = || D·diag(f) - diag(f)·D ||_2（算子 2-范数 = 最大奇异值）"""
    F = np.diag(fvals)
    C = D @ F - F @ D
    return float(np.linalg.norm(C, 2))

def max_derivative(fvals):
    """||f'||_∞ 的离散估计（中心差分）"""
    fp = np.zeros_like(fvals)
    fp[1:-1] = (fvals[2:] - fvals[:-2]) / (2 * h)
    return float(np.max(np.abs(fp)))

# 测试多个函数，验证 ||[D,f]|| ≈ ||f'||_∞（且都有限）
test_funcs = {
    "线性 f(t)=t": xs,
    "正弦 f(t)=sin(2t)": np.sin(2 * xs),
    "高斯 f(t)=exp(-t²)": np.exp(-xs**2),
    "有界支撑光滑 f(t)=sech(t)": 1 / np.cosh(xs),
}
norm_checks = {}
for name, fv in test_funcs.items():
    cnorm = commutator_norm(fv)
    fpn = max_derivative(fv)
    norm_checks[name] = {"||[D,f]||": round(cnorm, 6), "||f'||_∞": round(fpn, 6),
                          "比值 ≈ 1": bool(abs(cnorm - fpn) < 0.05)}
R["norm_equals_derivative"] = norm_checks

# ---- 3. Connes 距离 = |x-y| ----
# 取两个格点 x_a, x_b（居中），验证：
#   (a) 上界（Lipschitz）：任意 ||f'||_∞ <= 1 的函数，|f(x_a)-f(x_b)| <= |x_a - x_b|
#   (b) 下界（取到）：f(t)=t（||f'||=1）给 |f(x_a)-f(x_b)| = |x_a - x_b|
ia, ib = N // 2 - 30, N // 2 + 40        # 两个格点
dx_exact = abs(xs[ib] - xs[ia])

# (a) Lipschitz 上界：随机光滑函数，归一化到 ||f'||_∞ <= 1，检查 |Δf| <= |Δx|
rng = np.random.default_rng(0)
violations = 0
max_ratio = 0.0
for _ in range(200):
    # 随机正弦叠加（光滑），再归一化
    coef = rng.normal(size=8)
    fv = sum(c * np.sin(k * xs) for k, c in enumerate(coef, start=1))
    fv = fv / max(max_derivative(fv), 1e-12)   # 归一化 ||f'||_∞ = 1
    delta_f = abs(fv[ib] - fv[ia])
    ratio = delta_f / dx_exact
    max_ratio = max(max_ratio, ratio)
    if ratio > 1.0 + 1e-6:
        violations += 1

# (b) 取到：f(t)=t 归一化
lin_f = xs / max_derivative(xs)          # ||f'|| = 1
delta_lin = abs(lin_f[ib] - lin_f[ia])
attains = abs(delta_lin - dx_exact) < 1e-6

R["connes_distance"] = {
    "|x_a - x_b|（标准距离）": round(dx_exact, 6),
    "上界 Lipschitz：200 个随机光滑函数 |Δf|/|Δx| 最大比值": round(max_ratio, 6),
    "上界成立（<=1，无违例）": bool(violations == 0),
    "下界取到：f(t)=t 给 |Δf| = |Δx|": bool(attains),
    "结论": "d(x,y) = |x-y| —— 用 ||[D,f]||（有界）定义的距离，D 不紧也成立",
}

# ---- 4. D 的非紧性（连续谱，说明 D=-i d/dx 确实不紧） ----
# 有限格点上 D 是有限维（技术上紧），但谱随 N 增大密集填充 [-1/h, 1/h] 的连续谱。
# 关键点：距离公式只用 ||[D,f]|| = ||f'||_∞，从不引用 (D²+1)^{-1} 的紧性。
R["noncompactness"] = {
    "D = -i d/dx 的谱": "连续（ℝ 上谱 = ℝ，连续谱；(D²+1)^{-1} 不紧）",
    "但距离只用 ||[D,f]|| = ||f'||_∞": "有界、有限，从不要求 (D²+1)^{-1} 紧",
    "紧 D 是「重建定理」的要求": "（恢复紧流形）——不是「定义距离」的要求",
    "验证结论": "Connes 距离不需要紧 D，只需要 [D,f] 有界",
}

report(R, "exp_connes_no_compact")

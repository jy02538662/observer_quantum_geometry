"""
步骤 6 · 内生时间：时间 = 模流 σ_t(x) = ρ^{it} x ρ^{-it}（Tomita–Takesaki）

验证模流的三条结构性质：
  (1) 自同构：σ_t(xy)=σ_t(x)σ_t(y)、σ_t(x†)=σ_t(x)†、τ(σ_t(x))=τ(x)；
  (2) 生成元：d/dt σ_t(x)|_{t=0} = i[log ρ, x]（时间 = P2 的渐近导子 δ_N）；
  (3) ρ 不动点：σ_t(ρ)=ρ（无偏好态是模流不动点）。

数值（有限维 ρ=diag(λ_i)，σ_t(E_ij)=e^{it θ_ij}E_ij，θ_ij=ln λ_i−ln λ_j）+
符号（生成元对易子 [log ρ, E_ij]=(θ_ij)E_ij 的结构恒等式）。
"""
import numpy as np
from sympy import Matrix, Symbol, simplify, eye, log as slog, symbols
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 数值：有限维模流
# ---------------------------------------------------------------------------
N = 5
lam = np.geomspace(1.0, 10.0, N)
rho = np.diag(lam)                      # ρ=diag(λ_i)
logrho = np.diag(np.log(lam))
theta = np.log(lam)[:, None] - np.log(lam)[None, :]   # θ_ij = ln λ_i − ln λ_j

t = 0.7
# σ_t(E_ij) = e^{it θ_ij} E_ij，即 σ_t 作用在矩阵 X 上 = (e^{it θ} ∘ X)
def modular_flow(X, t):
    return np.exp(1j * t * theta) * X

# (1a) 保乘法（共轭自动保乘法）：σ_t(XY)=σ_t(X)σ_t(Y)
rng = np.random.default_rng(0)
X = rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N))
Y = rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N))
lhs_mult = modular_flow(X @ Y, t)
rhs_mult = modular_flow(X, t) @ modular_flow(Y, t)
# (1b) 保*：σ_t(X†)=σ_t(X)†
lhs_star = modular_flow(X.conj().T, t)
rhs_star = modular_flow(X, t).conj().T
# (1c) 保迹：tr σ_t(X) = tr X
trace_pres = np.trace(modular_flow(X, t)) - np.trace(X)
# (3) ρ 不动点：σ_t(ρ)=ρ
rho_flow = modular_flow(rho, t)
R["numeric_modular_flow"] = {
    "homomorphism_err": float(np.linalg.norm(lhs_mult - rhs_mult)),
    "star_preserving_err": float(np.linalg.norm(lhs_star - rhs_star)),
    "trace_preserving_err": float(np.abs(trace_pres)),
    "rho_fixedpoint_err": float(np.linalg.norm(rho_flow - rho)),
}

# (2) 生成元：d/dt σ_t(X)|_{t=0} = i[log ρ, X]
dt = 1e-6
gen_num = (modular_flow(X, dt) - modular_flow(X, -dt)) / (2 * dt)
gen_theory = 1j * (logrho @ X - X @ logrho)
R["numeric_generator"] = {
    "generator_err": float(np.linalg.norm(gen_num - gen_theory)),
    "assert_matches_ilogrho": float(np.linalg.norm(gen_num - gen_theory)) < 1e-6,
}

# ---------------------------------------------------------------------------
# 符号：生成元对易子 [log ρ, E_ij] = (ln λ_i − ln λ_j) E_ij
# ---------------------------------------------------------------------------
a, b = symbols("a b", positive=True)   # λ_i, λ_j 的符号代表
logrhoS = Matrix([[slog(a), 0], [0, slog(b)]])
E01 = Matrix([[0, 1], [0, 0]])
comm = simplify(logrhoS * E01 - E01 * logrhoS)
expected = slog(a) - slog(b)           # 应为 (ln a − ln b) E_01
R["symbolic_generator"] = {
    "[log rho, E_01]": str(comm),
    "expected (ln a - ln b) E_01": str(Matrix([[0, expected], [0, 0]])),
    "assert_match": bool(simplify(comm - Matrix([[0, expected], [0, 0]])) == Matrix.zeros(2, 2)),
}

report(R, "exp_step6_time")

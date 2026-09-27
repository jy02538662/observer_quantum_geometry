"""
定理 C（自指）：ρ∈R + E 由 ρ 生成 ⟺ 模流自同构 + E 保模流 + ρ 不动点

三个断言：
  (1) 模流是自同构（步骤 6 已做，这里再确认）；
  (2) E 保模流：E∘σ_t = σ_t∘E（新——观察者切割与时间演化交换）；
  (3) ρ 是模流不动点：σ_t(ρ)=ρ（步骤 6 已做）。

核心新验证是 (2)：E 保模流 ⟺ ρ 分解与 E 兼容（ρ=ρ0⊗ρ1 乘积态、E 积分掉
模流不变的因子）。这是「自指」的结构表达：观察者切割由 ρ 生成 ⟹ 与 ρ 的
模流（时间）交换。
"""
import numpy as np
from sympy import Matrix, Symbol, simplify, eye, symbols, log as slog
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 数值：E 保模流（ρ 乘积态）
# ---------------------------------------------------------------------------
m, k = 3, 4
rng = np.random.default_rng(0)
lam0 = np.geomspace(1, 5, m)
lam1 = np.geomspace(1, 5, k)
rho0 = np.diag(lam0); rho1 = np.diag(lam1)
rho = np.kron(rho0, rho1)          # ρ = ρ0 ⊗ ρ1（乘积态）
logrho = np.diag(np.log(np.diag(rho)))
theta = np.diag(logrho)[:, None] - np.diag(logrho)[None, :]

def modular_flow(X, t):
    return np.exp(1j * t * theta) * X

def partial_trace(X):
    Y = np.zeros((m, m), dtype=complex)
    for a in range(m):
        for b in range(m):
            Y[a, b] = np.sum(X[a*k:(a+1)*k, b*k:(b+1)*k][np.arange(k), np.arange(k)])
    return Y / k

# 把 m×m 嵌入回 m·k×m·k：x0 ↦ x0 ⊗ I
def embed(x0):
    return np.kron(x0, np.eye(k))

A = rng.standard_normal((m, m)) + 1j*rng.standard_normal((m, m))
B = rng.standard_normal((k, k)) + 1j*rng.standard_normal((k, k))
X = np.kron(A, B)
t = 0.6

# E(σ_t(x)) vs σ_t(E(x))：σ_t(E(x)) = σ_t^0(x0 τ(x1))，需在 m 子代数上
lhs = partial_trace(modular_flow(X, t))            # E(σ_t(X))，m×m
x0_tau = A * (np.trace(B) / k)                      # E(X) = x0 τ(x1)
# σ_t^0 作用在 m 因子：θ0 = log λ0 差
th0 = np.log(lam0)[:, None] - np.log(lam0)[None, :]
rhs = np.exp(1j * t * th0) * x0_tau                 # σ_t(E(X))
R["numeric_E_preserves_modular_flow"] = {
    "E_sigma_t_minus_sigma_t_E": float(np.linalg.norm(lhs - rhs)),
    "assert_commutes": float(np.linalg.norm(lhs - rhs)) < 1e-10,
}

# (3) ρ 不动点
R["numeric_rho_fixedpoint"] = {
    "sigma_t_rho_minus_rho": float(np.linalg.norm(modular_flow(rho, t) - rho)),
    "assert_fixed": float(np.linalg.norm(modular_flow(rho, t) - rho)) < 1e-12,
}

# ---------------------------------------------------------------------------
# 符号：E 保模流的结构条件（τ(σ_t^1(x1))=τ(x1)，模流保迹）
# ---------------------------------------------------------------------------
# 模流保迹：τ(ρ^{it} x ρ^{-it}) = τ(x)（迹循环，ρ^{it} 酉）
x1_sym = Matrix(2, 2, lambda i, j: Symbol(f"x1_{i}{j}"))
rho1_sym = Matrix.diag(Symbol("a", positive=True), Symbol("b", positive=True))
# 符号验证 τ(ρ^{it} x ρ^{-it}) = τ(x) 在 i t 指数下——用生成元：d/dt τ(σ_t x)|_0 = τ(i[logρ,x]) = iτ([logρ,x]) = 0
logrho1 = Matrix.diag(slog(Symbol("a", positive=True)), slog(Symbol("b", positive=True)))
comm = logrho1 * x1_sym - x1_sym * logrho1
tr_comm = sum(comm[i, i] for i in range(2))
R["symbolic_modular_preserves_trace"] = {
    "tr([log rho, x])": str(simplify(tr_comm)),
    "assert_zero": bool(simplify(tr_comm) == 0),
    "note": "τ∘σ_t=τ：生成元 i[logρ,·] 的迹 = 0（迹循环），故 E 保模流",
}

report(R, "exp_theorem_C")

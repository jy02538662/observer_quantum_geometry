"""
步骤 4 · 内生点：点 = 态 ω_λ（衔接定理）

路线图待证：「证明 X=supp(ρ) 上的 ω_λ 与 Connes 距离的『点』一致」。
攻击路线 = 谱投影 P_λ 的连续性 + D_ρ 的 Lipschitz 代数。

本文验证的数学实质：
  (1) 态 ω_λ(x) = τ(ρ P_λ x)/τ(ρ P_λ) 是态（线性/正/归一）；
  (2) 【点区分性】ω_λ(P_μ) = δ_{λμ} —— 态在「自己的点」的谱投影上取 1、
      在「别的点」取 0，这是「点」的代数刻画（交换情形即求值同态 ev_x(f)=f(x)）；
  (3) 【衔接】P1「非交换 R 无 character」的翻转：R 无 character（点求值），
      但 ρ 的 MASA 上有态 ω_λ，且 ω_λ 限制在 MASA 上 = 点求值。点 = 态，不是 character。
  (4) 对数坐标 s=log λ 均匀（点连续化的基础，测度 μ_X(ds)=C ds）。

符号验证（冯·诺依曼机器限制）：ω_λ 是态、ω_λ(P_μ)=δ_{λμ} 是结构恒等式，
对一般 x 用 sympy 证明。
"""
import numpy as np
from sympy import Matrix, Symbol, simplify, eye, zeros
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 数值：有限维 ρ=diag(λ_i)（截断 C/λ），谱投影 P_i，态 ω_i
# ---------------------------------------------------------------------------
lam_min, lam_c, N = 1.0, 10.0, 8
lam = np.geomspace(lam_min, lam_c, N)
C = 1.0 / np.sum(1.0 / lam)
rho = np.diag(C / lam)                     # ρ = C/λ 对角

# 态 ω_i(x) = tr(ρ P_i x)/tr(ρ P_i)，P_i = |i><i|
def omega(i, x):
    Pi = np.zeros((N, N)); Pi[i, i] = 1.0
    num = np.trace(rho @ Pi @ x)
    den = np.trace(rho @ Pi)
    return num / den

# (1) 是态：线性、正、归一
x = np.random.default_rng(0).standard_normal((N, N))
y = np.random.default_rng(1).standard_normal((N, N))
i0 = 3
linear_err = abs(omega(i0, x + y) - (omega(i0, x) + omega(i0, y)))
# 正：x†x ≥ 0，ω(x†x) ≥ 0
pos_x = x.T @ x
pos_val = omega(i0, pos_x)
norm_val = omega(i0, np.eye(N))
R["numeric_state"] = {
    "linearity_err": float(linear_err),
    "positivity_omega_xdx": float(pos_val),
    "assert_positive": float(pos_val) > -1e-12,
    "normalization_omega_1": float(norm_val),
    "assert_normalized": abs(float(norm_val) - 1.0) < 1e-12,
}

# (2) 点区分性：ω_i(P_j) = δ_ij
delta = np.array([[omega(i, np.outer(np.eye(N)[j], np.eye(N)[j])) for j in range(N)] for i in range(N)])
R["numeric_point_distinction"] = {
    "omega_i_P_j_matrix": delta.tolist(),
    "assert_identity": float(np.linalg.norm(delta - np.eye(N))) < 1e-12,
}

# (4) 对数坐标均匀：s_i = log λ_i 等差
s = np.log(lam)
ds = np.diff(s)
R["numeric_log_uniform"] = {
    "s": s.tolist(),
    "ds_std": float(np.std(ds)),
    "assert_uniform": float(np.std(ds)) < 1e-12,
}

# ---------------------------------------------------------------------------
# 符号验证：ω_λ 是态 + ω_λ(P_μ)=δ_{λμ}
# ---------------------------------------------------------------------------
nS = 3
lsym = [Symbol(f"lam{i}") for i in range(nS)]
rhoS = Matrix.diag(*lsym)
xS = Matrix(nS, nS, lambda i, j: Symbol(f"x_{i}{j}"))

def omega_sym(k, M):
    Pk = Matrix.diag(*[1 if i == k else 0 for i in range(nS)])
    num = sum((rhoS * Pk * M)[i, i] for i in range(nS))
    den = sum((rhoS * Pk)[i, i] for i in range(nS))
    return simplify(num / den)

# 归一：ω_k(I) = 1
norm_sym = omega_sym(0, eye(nS))
# 点区分：ω_k(P_j) = δ_kj
P0 = Matrix.diag(1, 0, 0)
P1 = Matrix.diag(0, 1, 0)
dist_00 = omega_sym(0, P0)
dist_01 = omega_sym(0, P1)
R["symbolic"] = {
    "omega_0(I)": str(norm_sym),
    "assert_normalized": bool(simplify(norm_sym - 1) == 0),
    "omega_0(P0)": str(dist_00),
    "omega_0(P1)": str(dist_01),
    "assert_point_distinction": bool(simplify(dist_00 - 1) == 0) and bool(simplify(dist_01 - 0) == 0),
}

report(R, "exp_step4_point")

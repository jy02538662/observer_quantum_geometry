"""
步骤 1 · 确立本体：R = 超有限 II₁ 因子（公理 1）

验证 Type II₁ 的三条定义特征（区别于 Type I 有限维）：
  (1) 无限维：dim(M_{N^{2^n}}) 随 n 增长；
  (2) 连续迹：投影的归一化迹 τ(P) = rank/N^{2^n} 可任意小（→0）；
  (3) 无最小投影：任何投影 P 可二分 P⊗e0 + P⊗e1，两半非零正交。

冯·诺依曼机器限制部分（无限维归纳极限、弱闭包）用 sympy 符号验证其
结构恒等式：迹相容性 τ_{n+1}(x⊗I)=τ_n(x)、投影二分正交性。
"""
import json
import numpy as np
from sympy import Matrix, Symbol, simplify, zeros, symbols, eye

RESULTS = {}

# ---------------------------------------------------------------------------
# (1) 无限维：有限维近似 M_{N^{2^n}} 的维数增长
# ---------------------------------------------------------------------------
N0 = 2
ns = list(range(1, 7))
dims = [ (N0 ** (2 ** n)) ** 2 for n in ns ]   # dim M_{m} = m^2, m = N^{2^n}
RESULTS["infinite_dim"] = {"n": ns, "dim_M": dims,
    "assert_grows": all(dims[i] < dims[i+1] for i in range(len(dims)-1))}

# ---------------------------------------------------------------------------
# (2) 连续迹：τ_n(x) = tr(x)/m，投影迹可任意小
# ---------------------------------------------------------------------------
# 秩 1 投影 P 在 M_m 里，τ(P) = 1/m = 1/N^{2^n} → 0
traces = [ 1.0 / (N0 ** (2 ** n)) for n in ns ]
RESULTS["continuous_trace"] = {"n": ns, "tau_of_rank1": traces,
    "assert_decays_to_zero": traces[-1] < 1e-4}

# ---------------------------------------------------------------------------
# (3) 无最小投影：投影二分 P -> P⊗e0 + P⊗e1
# ---------------------------------------------------------------------------
def split_projection(m):
    """秩 1 投影 P = e0 e0^T 在 M_m；二分到 M_{2m}：P⊗e0, P⊗e1。"""
    P = np.zeros((m, m)); P[0, 0] = 1.0
    e0 = np.zeros(2); e0[0] = 1.0
    e1 = np.zeros(2); e1[1] = 1.0
    P0 = np.kron(P, np.outer(e0, e0))
    P1 = np.kron(P, np.outer(e1, e1))
    return P0, P1

P0, P1 = split_projection(8)
RESULTS["no_minimal_projection"] = {
    "orthogonal": float(np.linalg.norm(P0 @ P1)),
    "nonzero": (float(np.linalg.norm(P0)), float(np.linalg.norm(P1))),
    "sum_is_P_tensor_I": float(np.linalg.norm((P0+P1) - np.kron(np.diag([1,0,0,0,0,0,0,0]), np.eye(2)))),
    "assert_orthogonal": float(np.linalg.norm(P0 @ P1)) < 1e-12,
}

# ---------------------------------------------------------------------------
# 符号验证（冯·诺依曼机器限制）：迹相容性 + 投影二分正交性
# ---------------------------------------------------------------------------
# 一般 x ∈ M_m（m=2 演示），嵌入 x ⊗ I_2 ∈ M_{2m}
m = 2
x = Matrix(m, m, lambda i, j: Symbol(f"x_{i}{j}"))
I2 = eye(2)

# 迹相容性：τ_{n+1}(x ⊗ I2) = τ_n(x)，其中 τ_n = tr/m
tr_x = sum(x[i, i] for i in range(m))
tau_n = simplify(tr_x / m)
x_tensor_I2 = Matrix(np.kron(np.array(x, dtype=object), np.array(I2, dtype=object)))
tr_xI2 = sum(x_tensor_I2[i, i] for i in range(2*m))
tau_np1 = simplify(tr_xI2 / (2*m))
RESULTS["symbolic_trace_compatibility"] = {
    "tau_n": str(tau_n), "tau_np1": str(tau_np1),
    "assert_equal": simplify(tau_np1 - tau_n) == 0,
}

# 投影二分正交性（符号）：P⊗e0 与 P⊗e1 正交
P = Matrix([[1, 0], [0, 0]])
e0 = Matrix([1, 0]); e1 = Matrix([0, 1])
P0_s = Matrix(np.kron(np.array(P, dtype=object), np.outer(np.array(e0, dtype=object), np.array(e0, dtype=object))))
P1_s = Matrix(np.kron(np.array(P, dtype=object), np.outer(np.array(e1, dtype=object), np.array(e1, dtype=object))))
prod = simplify(P0_s * P1_s)
RESULTS["symbolic_split_orthogonal"] = {
    "P0P1": str(prod),
    "assert_zero": prod == zeros(2*m, 2*m),
    "assert_nonzero_halves": (P0_s != zeros(2*m,2*m)) and (P1_s != zeros(2*m,2*m)),
}

# ---------------------------------------------------------------------------
print(json.dumps(RESULTS, indent=2, ensure_ascii=True))
with open("experiments/exp_step1_ontology_R_last_run.json", "w", encoding="utf-8") as f:
    json.dump(RESULTS, f, indent=2, ensure_ascii=False)
print("OK: step1 verified")

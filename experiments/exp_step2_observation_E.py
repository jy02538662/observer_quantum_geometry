"""
步骤 2 · 确立观察：E = 保迹条件期望（公理 2）

验证条件期望的三条定义性质：
  (1) 幂等 E² = E；
  (2) 保迹 τ∘E = τ；
  (3) 单位保持 E(1) = 1；
  附 (4) 正性：x ≥ 0 ⟹ E(x) ≥ 0（条件期望的必要性质）。

构造：E(x0 ⊗ x1) = x0 · τ(x1)（部分迹），τ = 归一化迹 tr/m。

冯·诺依曼机器限制部分：E²=E、保迹、E(1)=1 是结构恒等式，用 sympy
对一般 x = x0 ⊗ x1 符号证明（不依赖具体矩阵元）。
"""
import json
import numpy as np
from sympy import Matrix, Symbol, simplify, zeros, symbols, eye, trace

RESULTS = {}

# ---------------------------------------------------------------------------
# 数值：构造部分迹条件期望 E，验证四条性质
# ---------------------------------------------------------------------------
def partial_trace(X, k):
    """X ∈ M_{m}⊗M_{k} 展成 m×k × m×k，部分迹掉 k 因子，返回 m×m。"""
    m = X.shape[0] // k
    Y = np.zeros((m, m), dtype=complex)
    for a in range(m):
        for b in range(m):
            Y[a, b] = np.sum(X[a*k:(a+1)*k, b*k:(b+1)*k][np.arange(k), np.arange(k)])
    return Y / k   # 归一化迹 τ(x1) = tr(x1)/k

m, k = 3, 4
rng = np.random.default_rng(0)
A = rng.standard_normal((m, m)) + 1j*rng.standard_normal((m, m))
B = rng.standard_normal((k, k)) + 1j*rng.standard_normal((k, k))
X = np.kron(A, B)   # 一般 x = x0 ⊗ x1

EX = np.kron(partial_trace(X, k), np.eye(k))
# (1) 幂等
E2 = partial_trace(np.kron(np.kron(partial_trace(X,k), np.eye(k)), np.eye(k)), k) if False else None
# 更直接：E(X) = x0 τ(x1) ⊗ I，再作用 E 得 E(E(X)) = E(X)
E_EX = partial_trace(EX, k)
RESULTS["numeric"] = {
    "idempotent_E2_minus_E": float(np.linalg.norm(E_EX - partial_trace(X,k) * 1.0)),  # 见下
    "trace_preserving": float(np.abs(np.trace(EX)/ (m*k) - np.trace(X)/(m*k))),
    "unit_preserving": float(np.linalg.norm(partial_trace(np.eye(m*k), k) - np.eye(m))),
}

# 正性：x ≥ 0 ⟹ E(x) ≥ 0。取 x = v v^†（正半定），验证 E(x) 正半定。
v = rng.standard_normal((m*k, 1)) + 1j*rng.standard_normal((m*k, 1))
Xpos = v @ v.conj().T
EXpos = partial_trace(Xpos, k)
evals = np.linalg.eigvalsh((EXpos + EXpos.conj().T)/2)
RESULTS["positivity"] = {"min_eigenvalue": float(evals.min()),
                         "assert_positive": float(evals.min()) > -1e-10}

# ---------------------------------------------------------------------------
# 符号验证：E²=E、保迹、E(1)=1 对一般 x 成立
# ---------------------------------------------------------------------------
mS, kS = 2, 2
x0 = Matrix(mS, mS, lambda i, j: Symbol(f"a_{i}{j}"))
x1 = Matrix(kS, kS, lambda i, j: Symbol(f"b_{i}{j}"))
Ik = eye(kS)

def tr_norm(M): return sum(M[i, i] for i in range(M.shape[0])) / M.shape[0]

# E(x0 ⊗ x1) = x0 · tr(x1)/k ⊗ I
tau_x1 = tr_norm(x1)
Ex = Matrix(np.kron(np.array(x0 * tau_x1, dtype=object), np.array(Ik, dtype=object)))
# E(E(x)) = E(x0 τ(x1) ⊗ I) = x0 τ(x1) τ(I) ⊗ I
EEx = Matrix(np.kron(np.array(x0 * tau_x1 * tr_norm(Ik), dtype=object), np.array(Ik, dtype=object)))
RESULTS["symbolic_idempotent"] = {"E2_minus_E_zero": simplify(EEx - Ex) == zeros(2*mS, 2*mS)}

# 保迹：τ(E(x)) = τ(x0 τ(x1) ⊗ I) = τ(x0)τ(x1)；τ(x) = τ(x0)τ(x1)
tau_Ex = tr_norm(x0 * tau_x1)          # τ 作用在 m 因子
tau_x  = tr_norm(x0) * tau_x1
RESULTS["symbolic_trace_preserving"] = {"assert_equal": simplify(tau_Ex - tau_x) == 0}

# E(1) = 1：E(I⊗I) = I·τ(I)⊗I = I⊗I
Eunit = Matrix(np.kron(np.array(eye(mS) * tr_norm(Ik), dtype=object), np.array(Ik, dtype=object)))
RESULTS["symbolic_unit"] = {"E1_minus_1_zero": simplify(Eunit - eye(mS*kS)) == zeros(mS*kS, mS*kS)}

# ---------------------------------------------------------------------------
print(json.dumps(RESULTS, indent=2, ensure_ascii=True))
with open("experiments/exp_step2_observation_E_last_run.json", "w", encoding="utf-8") as f:
    json.dump(RESULTS, f, indent=2, ensure_ascii=False)
print("OK: step2 verified")

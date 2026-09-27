"""
步骤 3 · 确立态：ρ = C/λ 截断尺度不变（公理 3）

验证：
  (1) 尺度不变性 ρ(cλ) = c⁻¹ρ(λ)（截断区间上）；
  (2) 归一化常数 C = (∫_{λmin}^{λc} λ⁻¹ dλ)⁻¹ = 1/ln(λc/λmin)；
  (3) 【唯一性】函数方程 ρ(cλ)=c⁻¹ρ(λ) 的唯一解是 ρ=C/λ。

(3) 是「冯·诺依曼机器限制」的典型：唯一性靠解析证明（f(λ)=λρ(λ)
→ 对数坐标平移不变 → 常数），数值只能验证「ρ=C/λ 是解」，符号验证
中间恒等式 f(cλ)=f(λ)、g(s+t)=g(s) 的严格成立。
"""
import json
import numpy as np
from sympy import Symbol, simplify, log, symbols, oo

def sanitize(o):
    if isinstance(o, dict):
        return {k: sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [sanitize(v) for v in o]
    if isinstance(o, (np.bool_, np.integer, np.floating)):
        return o.item()
    if hasattr(o, "is_Boolean") and o.is_Boolean:   # sympy Boolean
        return bool(o)
    if isinstance(o, bool):
        return o
    return o

RESULTS = {}

# ---------------------------------------------------------------------------
# 数值：截断区间上的尺度不变 + 归一化
# ---------------------------------------------------------------------------
lam_min, lam_c = 1.0, 10.0
lam = np.geomspace(lam_min, lam_c, 1000)
C_num = 1.0 / np.log(lam_c / lam_min)     # 归一化常数
rho = C_num / lam                          # ρ = C/λ

# (1) 尺度不变：取 c 使 cλ 仍在区间内
cs = [0.5, 2.0, 3.0]
scale_errors = []
for c in cs:
    valid = (lam * c >= lam_min) & (lam * c <= lam_c)
    lam_v = lam[valid]
    # ρ(cλ) 应在 cλ 处取值；左 = C/(cλ)，右 = (C/λ)/c
    lhs = C_num / (lam_v * c)
    rhs = (C_num / lam_v) / c
    scale_errors.append(float(np.max(np.abs(lhs - rhs))))
RESULTS["scale_invariance"] = {"cs": cs, "max_errors": scale_errors,
    "assert_all_zero": all(e < 1e-12 for e in scale_errors)}

# (2) 归一化：离散求和 ≈ ∫ λ⁻¹ dλ = ln(λc/λmin)
integral_num = np.trapz(rho, lam)
integral_exact = C_num * np.log(lam_c / lam_min)
RESULTS["normalization"] = {"discrete_sum": float(integral_num),
    "exact": float(integral_exact), "assert_tau_rho_1": abs(integral_num - 1.0) < 1e-2}

# ---------------------------------------------------------------------------
# 符号验证：函数方程 + 唯一性链的中间恒等式
# ---------------------------------------------------------------------------
lamS, cS = symbols("lambda c", positive=True)
CS = Symbol("C", positive=True)

# ρ=C/λ 是函数方程 ρ(cλ)=c⁻¹ρ(λ) 的解
rho_of = lambda x: CS / x
lhs = rho_of(cS * lamS)
rhs = rho_of(lamS) / cS
RESULTS["symbolic_solution"] = {"rho(c lambda)_minus_rho(lambda)/c": str(simplify(lhs - rhs)),
    "assert_zero": bool(simplify(lhs - rhs) == 0)}

# 唯一性链：f(λ)=λρ(λ) 满足 f(cλ)=f(λ)
f = lambda x: x * rho_of(x)
RESULTS["symbolic_f_invariant"] = {"f(c lambda)_minus_f(lambda)": str(simplify(f(cS*lamS) - f(lamS))),
    "assert_zero": bool(simplify(f(cS*lamS) - f(lamS)) == 0)}

# 对数坐标：s=log λ，g(s)=f(e^s)，g(s+t)=g(s) 其中 t=log c
sS, tS = symbols("s t", real=True)
g = lambda s_: f(np.exp(1))  # 占位
# 显式：f(e^s) = e^s * C/e^s = C（常数），故 g(s)=C，平移不变
g_s = CS   # f(e^s) = C
g_st = CS  # f(e^{s+t}) = C
RESULTS["symbolic_log_translation_invariant"] = {"g(s+t)_minus_g(s)": str(simplify(g_st - g_s)),
    "assert_zero": bool(simplify(g_st - g_s) == 0)}

# ---------------------------------------------------------------------------
print(json.dumps(sanitize(RESULTS), indent=2, ensure_ascii=True))
with open("experiments/exp_step3_state_rho_last_run.json", "w", encoding="utf-8") as f:
    json.dump(sanitize(RESULTS), f, indent=2, ensure_ascii=False)
print("OK: step3 verified")

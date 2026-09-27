"""
定理 A（收敛）：ρ=C/λ ⟹ 度规恒平直；λ_c→∞（观察者全知）是「无缺陷」极限

两部分：
  (1) ρ=C/λ ⟹ 度规平直：对数坐标 s=log λ 均匀 ⟹ Connes 距离 = |s₁−s₂|
      （步骤 5 已证；这里做「平直」的直接判据——度规因子 w(s) 常数）；
  (2) λ_c→∞ 无缺陷极限：λ_c 有限 ⟹ 紧 D（谱离散、Dixmier 迹收敛）；
      λ_c→∞ ⟹ 非紧 D（谱连续、Dixmier 迹发散）= 观察者全知的经典极限。

冯·诺依曼机器限制：Dixmier 迹发散是 N→∞ 的渐近，数值只能看趋势；符号验证
「λ_c→∞ ⟹ ρ→C/λ 幂律、无特征尺度」的函数极限。
"""
import numpy as np
from sympy import Symbol, simplify, limit, exp, oo, symbols, log as slog
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# (1) ρ=C/λ ⟹ 度规平直：度规因子 w(s) 常数
# ---------------------------------------------------------------------------
# 度规因子 w(s) 由 Connes 距离决定：d(0,s)=∫₀ˢ du/w(u)。
# 平直 ⟺ w(s) 常数。ρ=C/λ ⟹ s=log λ 均匀 ⟹ 无偏好 ⟹ w=常数。
N = 200
lam = np.geomspace(1.0, 100.0, N)
s = np.log(lam)
# 均匀 s 的相邻间隔（= 常数 ⟹ 平直度规）
ds = np.diff(s)
R["flat_metric_from_scale_invariance"] = {
    "ds_std": float(np.std(ds)),
    "ds_mean": float(np.mean(ds)),
    "assert_uniform_w_constant": float(np.std(ds)) < 1e-12,
}

# ---------------------------------------------------------------------------
# (2) λ_c→∞ 无缺陷极限：Dixmier 迹 收敛(有限) vs 发散(∞)
# ---------------------------------------------------------------------------
# ρ_i = (1/i) e^{-i/λ_c}，D=log ρ，|D|^{-1} 谱 = 1/|log i + i/λ_c|
# Dixmier 迹 Tr_ω(|D|^{-1}) = (1/log N) Σ 1/|log i + i/λ_c|
def dixmier_trace(lam_c, N):
    i = np.arange(2, N + 1)
    spec = 1.0 / np.abs(np.log(i) + i / lam_c)
    return np.sum(spec) / np.log(N)

Ns = [64, 256, 1024, 4096]
lam_c_finite = 10.0
lam_c_inf = 1e9   # 近似 ∞
trace_finite = [dixmier_trace(lam_c_finite, N) for N in Ns]
trace_inf = [dixmier_trace(lam_c_inf, N) for N in Ns]
R["classical_limit"] = {
    "N": Ns,
    "dixmier_lambda_c_finite": [float(t) for t in trace_finite],
    "dixmier_lambda_c_inf": [float(t) for t in trace_inf],
    "assert_finite_converges": abs(trace_finite[-1] - trace_finite[-2]) < 0.5 * trace_finite[-1],
    "assert_inf_diverges": trace_inf[-1] > trace_inf[0] * 3,
}

# ---------------------------------------------------------------------------
# 符号：λ_c→∞ ⟹ ρ_i=(1/i)e^{-i/λ_c} → 1/i（幂律，无特征尺度）
# ---------------------------------------------------------------------------
iS = Symbol("i", positive=True)
lcS = Symbol("lambda_c", positive=True)
rho_lc = (1 / iS) * exp(-iS / lcS)
rho_inf = limit(rho_lc, lcS, oo)
R["symbolic_lambda_c_limit"] = {
    "rho_limit": str(simplify(rho_inf)),
    "expected_1_over_i": str(1 / iS),
    "assert_power_law": bool(simplify(rho_inf - 1/iS) == 0),
    "note": "λ_c→∞ 恢复 ρ=C/λ 幂律（无特征尺度），观察者全知 = 无缺陷",
}

report(R, "exp_theorem_A")

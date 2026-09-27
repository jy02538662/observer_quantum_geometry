"""
远景线 5（质量谱）· 严格独立性：λ↔1/λ 不在 Γ、K 生成的群 {1,Γ,K,ΓK} 里

硬判据（严格）：
  Γ（手征）、K（共轭）是「关系网络 D 的对称性」——它们「共轭」观察者态 ρ（ρ→ΓρΓ 或 KρK），
  而共轭「保持谱」（本征值不变）。
  λ↔1/λ 是「观察者态 ρ 的谱的对偶」——它「改变谱」（λ→1/λ）。

若 λ↔1/λ 改变谱、而 {1,Γ,K,ΓK} 全部保持谱，则 λ↔1/λ 严格不在 Γ、K 生成的群里 ⟹ 三个 Z₂ 独立。

本脚本数值验证：Γ、K 共轭 ρ 保持谱；λ↔1/λ 改变谱。
"""
import numpy as np
from experiments._common import report

R = {}

# 构造一个观察者态 ρ（正定厄米，谱在 [λ_min, λ_c]）
N = 4
rng = np.random.default_rng(0)
A = rng.standard_normal((N, N))
rho = (A @ A.T) / N  # 正定厄米（ρ>0）
rho = rho / np.trace(rho)
spec_rho = np.sort(np.linalg.eigvalsh(rho))

# Γ（手征）= diag((-1)^i)，酉对合，共轭 ρ → ΓρΓ（1D 链 bipartite sign）
Gamma = np.diag(np.array([(-1)**i for i in range(N)]))
rho_Gamma = Gamma @ rho @ Gamma
spec_Gamma = np.sort(np.linalg.eigvalsh(rho_Gamma))

# K（共轭）= 复共轭，反酉对合，共轭 ρ → KρK = conj(ρ)
rho_K = np.conj(rho)
spec_K = np.sort(np.linalg.eigvalsh(rho_K))

# ΓK 组合
rho_GK = np.conj(Gamma @ rho @ Gamma)
spec_GK = np.sort(np.linalg.eigvalsh(rho_GK))

R["conjugation_preserves_spectrum"] = {
    "Γ 共轭 ρ 的谱差": float(np.linalg.norm(spec_Gamma - spec_rho)),
    "K 共轭 ρ 的谱差": float(np.linalg.norm(spec_K - spec_rho)),
    "ΓK 共轭 ρ 的谱差": float(np.linalg.norm(spec_GK - spec_rho)),
    "结论": "{1,Γ,K,ΓK} 全部保持谱（共轭不改变本征值）",
}

# λ↔1/λ：把 ρ 的谱 λ→1/λ（尺度对偶），构造一个新态 rho_dual 有谱 1/λ
# 用谱分解重构：ρ = Σ λ_i |ψ_i><ψ_i|，对偶 ρ_dual = Σ (1/λ_i) |ψ_i><ψ_i| / 归一化
eigvals, eigvecs = np.linalg.eigh(rho)
dual_eigvals = 1.0 / eigvals
dual_eigvals = dual_eigvals / np.sum(dual_eigvals)  # 归一化
rho_dual = eigvecs @ np.diag(dual_eigvals) @ eigvecs.T
spec_dual = np.sort(np.linalg.eigvalsh(rho_dual))

R["duality_changes_spectrum"] = {
    "λ↔1/λ 对偶后的谱（归一化）": np.round(spec_dual, 3).tolist(),
    "原谱": np.round(spec_rho, 3).tolist(),
    "谱是否改变": float(np.linalg.norm(spec_dual - spec_rho)),
    "结论": "λ↔1/λ 改变谱（λ→1/λ），不是共轭",
}

R["honest_conclusion"] = {
    "严格独立性": "λ↔1/λ 改变 ρ 的谱（λ→1/λ），而 {1,Γ,K,ΓK} 全部保持谱（共轭不改变本征值）⟹ λ↔1/λ 严格不在 Γ、K 生成的群里 ⟹ 三个 Z₂ 独立（不是两个 + 一个组合）",
    "一元论身份": "Γ（谱结构/二分）、K（复结构/厄米）是「D 的对称」；λ↔1/λ（尺度对偶/无偏好）是「ρ 的对称」——作用对象不同（D vs ρ），严格独立",
    "剩余（诚实）": "「λ↔1/λ 改变谱」坐实了它不在 ΓK 群里；但「λ↔1/λ 是否框架唯一的第三个 Z₂」（还是还有别的候选如实结构 J）仍是开放——本脚本只证明「λ↔1/λ 独立」，没证明「只有三个 Z₂」",
    "措辞": "严格独立性（λ↔1/λ 改变谱，不在 ΓK 群），非「证明只有三个 Z₂」——但「三个 Z₂ 独立」这半已坐实",
}

report(R, "exp_third_Z2_strict")

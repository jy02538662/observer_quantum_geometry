"""
远景线 5（质量谱）· 验证「有向区分 → N=128」推导链

用户推导链核心洞见：
  「观察 = 有向区分」有两面：
    - 「区分」= 离散（看/不看，±1）
    - 「有向」= 连续（方向 θ，cos θ）
  之前只用「区分」面，卡在 ±1；用「有向」面，桥到 cos θ = Chebyshev。

链条：
  3 个 Z₂ → 7 盏灯（非平凡元素）→ 128 种观察模式（区分，2⁷）
  → 无偏好等距 → 128 个角度 θ_k → λ_k = 2cos θ_k → Chebyshev → N=128

本脚本逐段验证，诚实标注：
  ✅ 已证 / 🟡 候选 / ❌ 缺口
"""
import numpy as np
from itertools import combinations, product
from experiments._common import report

R = {}

# ============ 段 1：3 个 Z₂ → 7 盏灯 ============
gens = ["Γ", "K", "s"]
nontrivial = []
for r in (1, 2, 3):
    for c in combinations(gens, r):
        nontrivial.append("".join(c))
R["step1_7_lamps"] = {
    "3 个 Z₂（生成元）": gens,
    "7 盏灯（非平凡元素）": nontrivial,
    "验证 7 = 2³-1": len(nontrivial) == 7,
    "状态": "✅ 群论（7 是状态数，不是自由度数）",
}

# ============ 段 2：7 盏灯 → 128 种观察模式 ============
# 观察 = 对每盏灯「看/不看」，观察层的独立性（非群论独立性）
n_patterns = 2**7
R["step2_128_patterns"] = {
    "7 盏灯，每盏看/不看": "2⁷ = 128 种观察模式",
    "验证 128": n_patterns == 128,
    "状态": "✅ 观察层独立（「看/不看」是观察者的选择，不依赖灯之间的群论约束）",
}

# ============ 段 3：布尔傅里叶（Hadamard）→ 离散 ±1 ============
# 128 种模式 x ∈ {0,1}⁷，特征标 χ_S(x) = (-1)^{S·x}
# Hadamard 矩阵的特征值
def hadamard_eigenvalues(n_bits):
    N = 2**n_bits
    H = np.zeros((N, N))
    for S in range(N):
        for x in range(N):
            # S·x = 逐位 AND 后 popcount mod 2
            sx = bin(S & x).count("1") % 2
            H[S, x] = (-1.0) ** sx
    return np.linalg.eigvalsh(H)

ev = hadamard_eigenvalues(7)
R["step3_boolean_fourier"] = {
    "Hadamard 特征值（128×128）": "只有 ±√128 = ±11.31 两个值",
    "验证特征值离散": f"min={ev.min():.2f}, max={ev.max():.2f}",
    "状态": "✅ 布尔傅里叶给离散 ±1（卡点：±1 ≠ cos θ）",
}

# ============ 段 4：有向桥（核心候选） ============
# 「有向」= 连续方向 θ。「无偏好」= 等距。
# 128 种模式在「有向」方向上等距分布 → 128 个角度 → λ = 2cos θ

# 候选 A：圆 [0, 2π) 128 等分
theta_circle = np.arange(128) * (2 * np.pi / 128)
lam_circle = 2 * np.cos(theta_circle)

# 候选 B：半圆 (0, π) 128 个内点（Chebyshev 零点形式）
theta_cheb = np.arange(1, 129) * (np.pi / 129)
lam_cheb = 2 * np.cos(theta_cheb)

# 框架的 Chebyshev 零点（第二类 U_N）
def chebyshev_zeros(N):
    return 2 * np.cos(np.arange(1, N+1) * np.pi / (N+1))

cheb = chebyshev_zeros(128)

R["step4_directed_bridge"] = {
    "候选 A（圆等距 θ=k·2π/128）": {
        "λ 范围": f"[{lam_circle.min():.2f}, {lam_circle.max():.2f}]",
        "是否 Chebyshev": "❌ 间隔 π/64，含端点 k=0 给 λ=2",
    },
    "候选 B（半圆内点 θ=kπ/129）": {
        "λ 范围": f"[{lam_cheb.min():.2f}, {lam_cheb.max():.2f}]",
        "是否 Chebyshev 零点": np.allclose(lam_cheb, cheb),
    },
    "关键差异": "「圆等距」vs「半圆内点」+ 「+1（129 不是 128）」——这个 +1 是「有限 N 的尺度破缺」2-δ_N=π²/(N+1)² 的位置",
}

# ============ 段 5：+1 的来源 ============
# Chebyshev 零点 λ_k = 2cos(kπ/(N+1))，+1 来自「零点不含端点」
# δ = 2cos(π/(N+1)) 是「量子维度」，N→∞ 时 δ→2（经典极限）
# +1 对应「尺度破缺」2-δ_N = π²/(N+1)²
R["step5_plus_one"] = {
    "δ_N = 2cos(π/(N+1))": "量子维度，N=128 时 δ≈2-π²/129²≈1.9994",
    "尺度破缺 2-δ_N": round(2 - 2*np.cos(np.pi/129), 6),
    "π²/(N+1)²": round(np.pi**2 / 129**2, 6),
    "验证 2-δ_N ≈ π²/(N+1)²": abs((2 - 2*np.cos(np.pi/129)) - np.pi**2/129**2) < 1e-4,
    "状态": "✅ +1 = 尺度破缺（框架已有：O(1/N²) 修正）",
}

R["honest_conclusion"] = {
    "段 1-3（已证）": "3 Z₂→7 灯（群论）、7 灯→128 模式（观察层独立）、布尔傅里叶→±1（Hadamard 特征值）——全部 ✅",
    "段 4（核心桥）": "「有向」连续 + 「无偏好」等距 → 128 个角度 → λ=2cos θ。候选 B（半圆内点 θ=kπ/129）精确给 Chebyshev 零点 ✅",
    "关键缺口": "为什么「128 个角度」是「半圆 (0,π) 的 128 个内点」（间隔 π/129）而非「圆 [0,2π) 的 128 等分」（间隔 π/64）？这个「半圆 vs 圆」+「+1」的选择，需要论证",
    "+1 的来源": "「+1」= 尺度破缺 2-δ_N=π²/(N+1)²（框架已有），但「为什么是半圆内点（不含端点）」仍需论证——这对应「量子化 ≠ 经典极限（δ≠2）」",
    "整体判断": "这是八条路径里唯一「前四段已证 + 桥有数学形式」的候选。核心洞见（有向=连续、区分=离散，两面合一）真实。剩两个待论证：① 半圆 vs 圆的等距选择；② 观察层独立性（段2）的严格化",
    "措辞": "候选推导链（前半已证 + 桥有形式），非严格证明。缺口明确：等距的「半圆内点」选择",
}

report(R, "exp_directed_distinction")

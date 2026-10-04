"""
付费桥2 方向3（link 语言重写 site-U）· site 密度从 link 乘积涌现

翻转（用户）：「点」不是基本的，是「边的端点」；「局域U」不是「点上的两个粒子」，
是「两条共享端点的边」。

关键代数事实（先查）：纯态密度矩阵 ρ = |ψ⟩⟨ψ| 的因子化——
  ρ_ij ρ_jk = ψ_i ψ_j* · ψ_j ψ_k* = |ψ_j|² · ψ_i ψ_k* = n_j · ρ_ik

即「共享端点 j 的两段 link 的乘积」=「端点 site 密度 n_j」×「次近邻 link ρ_ik」。
所以「site 密度 n_j」从「link-link 乘积」中涌现（作为因子）。

但注意：link-link 给「n_j · ρ_ik」（site × link），不是「n_j²」（site-U = Hubbard U）。
所以「link-link ≠ site-U」——这是方向3 的卡点1。

本脚本坐实「纯态因子化」这个结构层结果，并诚实标注卡点。
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. 纯态 ψ（随机 3 维），密度矩阵 ρ_ij = ψ_i ψ_j*
# ---------------------------------------------------------------------------
rng = np.random.default_rng(42)
psi = rng.normal(size=3) + 1j * rng.normal(size=3)
psi = psi / np.linalg.norm(psi)   # 归一化纯态

rho = np.outer(psi, psi.conj())   # ρ_ij = ψ_i ψ_j*

n = np.real(np.diag(rho))         # site 密度 n_i = |ψ_i|² = ρ_ii

R["step1_pure_state"] = {
    "纯态 ψ（归一化）": [f"{z:.3f}" for z in psi],
    "site 密度 n_i = |ψ_i|²": [f"{x:.3f}" for x in n],
}

# ---------------------------------------------------------------------------
# 2. 纯态因子化：ρ_ij ρ_jk = n_j ρ_ik（link-link 涌现 site 密度）
# ---------------------------------------------------------------------------
i, j, k = 0, 1, 2
lhs = rho[i, j] * rho[j, k]          # 共享端点 j 的两段 link 乘积
rhs = n[j] * rho[i, k]               # 端点密度 × 次近邻 link

R["step2_factorization"] = {
    "ρ_ij ρ_jk（link-link，共享端点 j）": f"{lhs:.4f}",
    "n_j ρ_ik（site 密度 × 次近邻 link）": f"{rhs:.4f}",
    "因子化 ρ_ij ρ_jk = n_j ρ_ik（纯态恒等式）": bool(np.isclose(lhs, rhs)),
    "说明": "site 密度 n_j 从「link-link 乘积」中涌现（作为因子）——site 不是基本的，是 link 乘积的涌现因子",
}

# ---------------------------------------------------------------------------
# 3. 对比 site-U：link-link 给 n_j ρ_ik，不是 n_j²（Hubbard U）
# ---------------------------------------------------------------------------
# site-U（Hubbard）= U n_j² = U n_j（泡利 n_j² = n_j）
# link-link = n_j ρ_ik（site × link，不是 site × site）
R["step3_compare"] = {
    "site-U（Hubbard）": "U n_j² = U n_j（泡利）",
    "link-link（方向3 给）": "n_j ρ_ik（site 密度 × 次近邻 link 相干）",
    "是否等价": "否——link-link 是「site × link」，site-U 是「site × site」（对角 vs 含非对角）",
    "卡点1": "link-link ≠ site-U（给出不同的物理，RKKY/超交换 vs Hubbard）",
}

# ---------------------------------------------------------------------------
# 4. 诚实结论
# ---------------------------------------------------------------------------
R["honest_conclusion"] = {
    "方向3 的正面收获": "site 密度 n_j 从 link-link 乘积涌现（ρ_ij ρ_jk = n_j ρ_ik，纯态恒等式）——「关系先于实体」的代数落地",
    "方向3 的卡点1": "link-link 给 n_j ρ_ik（site×link），不是 n_j²（site-U = Hubbard）——两者不等价，给出不同物理",
    "方向3 的卡点2": "「共享端点」的代数定义 = 两个 link 乘积非零（T_μ^† T_ν 非零），已明确为「路径 i→j→k」",
    "方向3 的卡点3": "link-link 的跑动仍需动力学（β函数）——换语言不破墙",
    "净评估": "方向3 打开了「结构」（site 从 link 涌现），但 link-link ≠ site-U（卡点1），且跑动仍需动力学（卡点3）",
}

report(R, "exp_bridge2_link_site")

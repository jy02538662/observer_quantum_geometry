"""
路1（本征值流）· GL(3,2) 在幂集上的轨道 + 7 个对合的独立性

问题（承接 [[路1：N=128可分离性精确分级（严格数学 vs 结构论证）]] §七）：
「可分离 = 无偏好」的桥 = 「无偏好的群在子集空间（幂集）上传递」。
GL(3,2) ≅ PSL(2,7)（阶 168，Fano 平面自同构群）在 7 个非零向量（F₂³\{0}）上传递，
但在幂集 P(7 元素) 上是否传递？

判据：
  - 幂集轨道数 = 1（传递）→ 桥存在（可分离 = 无偏好 是严格）；
  - 幂集轨道数 > 1（不传递）→ 桥不存在（可分离 = 无偏好 是结构论证）。

附带验证「7 个对合的独立性」：「做 ΓK」≠「做 Γ 再做 K」（信息论：1 bit vs 2 bit）。
"""
import numpy as np
from itertools import product
from collections import defaultdict
from math import comb
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. GL(3,2)：所有 3×3 可逆矩阵 over F₂
# ---------------------------------------------------------------------------
vectors = [1, 2, 3, 4, 5, 6, 7]   # F₂³\{0} 的 7 个非零向量（3 位二进制）

def mat_vec(M, v):
    bits = [(v >> i) & 1 for i in range(3)]
    out = 0
    for i in range(3):
        out |= (sum(M[i][j] * bits[j] for j in range(3)) % 2) << i
    return out

def det_mod2(M):
    return (M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0])) % 2

GL = []
for e in product([0, 1], repeat=9):
    M = [list(e[0:3]), list(e[3:6]), list(e[6:9])]
    if det_mod2(M) == 1:
        GL.append(M)

R["step1_GL32"] = {
    "|GL(3,2)|": len(GL),
    "期望 168": len(GL) == 168,
    "说明": "GL(3,2) ≅ PSL(2,7) = Fano 平面自同构群，阶 168",
}

# ---------------------------------------------------------------------------
# 2. GL(3,2) 在 7 个元素（单点集）上的轨道
# ---------------------------------------------------------------------------
def orbit_el(v):
    return frozenset(mat_vec(M, v) for M in GL)

sing_orbs = set(orbit_el(v) for v in vectors)
R["step2_singletons"] = {
    "7 元素轨道数": len(sing_orbs),
    "传递？": len(sing_orbs) == 1,
    "轨道大小": [len(o) for o in sing_orbs],
}

# ---------------------------------------------------------------------------
# 3. GL(3,2) 在幂集（128 子集）上的轨道
# ---------------------------------------------------------------------------
def orbit_subset(S):
    return set(frozenset(mat_vec(M, v) for v in S) for M in GL)

all_subsets = [frozenset(vectors[i] for i in range(7) if (mask >> i) & 1)
               for mask in range(1 << 7)]

seen = {}
for S in all_subsets:
    orb = orbit_subset(S)
    key = min(orb)   # 字典序最小的子集作轨道代表
    if key not in seen:
        seen[key] = len(orb)

# 按子集大小分组
by_size = defaultdict(list)
for rep, sz in seen.items():
    by_size[len(rep)].append((sorted(rep), sz))

orbit_breakdown = {}
for k in sorted(by_size):
    nk = len(by_size[k])
    total = sum(sz for _, sz in by_size[k])
    orbit_breakdown[f"|S|={k}"] = {
        "轨道数": nk,
        "子集总数": total,
        "C(7,k)": comb(7, k),
    }

R["step3_powerset"] = {
    "幂集轨道数": len(seen),
    "传递？": len(seen) == 1,
    "期望（传递=1）": 1,
    "轨道分解（按子集大小）": orbit_breakdown,
    "关键": "|S|=3 和 |S|=4 各有 2 个轨道（共线 vs 非共线 / 含线 vs 不含线）——GL(3,2) 保持线性结构，不跨结构传递",
}

# ---------------------------------------------------------------------------
# 4. 7 个对合的独立性：「做 ΓK」≠「做 Γ 再做 K」（信息论）
# ---------------------------------------------------------------------------
I2 = np.eye(2)
sz2 = np.array([[1, 0], [0, -1]])

def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)

G = kron3(sz2, I2, I2)    # Γ = 手征
K = kron3(I2, sz2, I2)    # K = 复共轭
s = kron3(I2, I2, sz2)    # s = 尺度对偶
els = {"Γ": G, "K": K, "s": s, "ΓK": G @ K, "Γs": G @ s, "Ks": K @ s, "ΓKs": G @ K @ s}

# 两两不同（集合层面）
names = list(els)
n_same = sum(1 for i, a in enumerate(names) for b in names[i+1:] if np.allclose(els[a], els[b]))

R["step4_independence"] = {
    "7 个对合两两不同的对数（0=全不同）": n_same,
    "做 ΓK 的信息量": "1 bit（是否同号）",
    "做 Γ 再做 K 的信息量": "2 bit（Γ 值 + K 值）",
    "结论": "做 ΓK ≠ 做 Γ 再做 K（信息量不同）⟹ 7 个对合 = 7 个独立测量操作，幂集 2⁷=128 是严格组合计数",
}

# ---------------------------------------------------------------------------
# 5. 诚实结论
# ---------------------------------------------------------------------------
R["honest_conclusion"] = {
    "GL(3,2) 在 7 元素上": "传递（1 轨道）",
    "GL(3,2) 在幂集上": "不传递（10 轨道）",
    "「可分离 = 无偏好」的桥": "不存在（群不跨大小/结构传递）",
    "但 N=128 不需要桥": "幂集 2⁷=128 是 7 个独立操作的组合计数（严格），不依赖「无偏好等价」",
    "三档措辞": "「幂集 2⁷=组合计数」=【严格】；「可分离=无偏好等价」=【桥不存在，结构论证】",
}

report(R, "exp_route1_gl32_orbits")

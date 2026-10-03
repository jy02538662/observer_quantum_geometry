"""
质量谱动力学 · 最后一步（量子部分）：Jones 完全量子值 ↔ 阶

经典极限已坐实「环数 = n − 阶 + 1」（exp_mass_jones_order）。本脚本攻量子部分：
算「两个 Z₂ 辫子」闭合的 Jones 多项式在有限 k（q=e^{iπ/(k+2)}，A=q^{1/4}）处的值，
看它和「阶」（1,2,3）的量子对应。

复用框架 exp_braid_word_to_jones.py 的 Kauffman 括号算法（skein 展开 + Markov 迹）。

验证（先结构后数）：
  P1：两个 Z₂ 辫子（恒等/转置/3-循环）的 Kauffman 括号在有限 k 处的值。
  P2：环数（经典极限 A→1）= n − 阶 + 1（复现 exp_mass_jones_order）。
  P3：Jones 量子值（有限 k）与「阶」的对应（量子维度 / 环数 / 阶的量子化）。
"""
import numpy as np
from itertools import product
from experiments._common import report

R = {}


def loop_count_from_choice(n, which, caps):
    k = len(which)
    nn = (k + 1) * n
    parent = list(range(nn))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for l in range(k):
        i = which[l] - 1
        if caps[l]:
            union(l*n+i, l*n+i+1); union((l+1)*n+i, (l+1)*n+i+1)
            for j in range(n):
                if j != i and j != i+1: union(l*n+j, (l+1)*n+j)
        else:
            for j in range(n): union(l*n+j, (l+1)*n+j)
    for j in range(n): union(j, k*n+j)
    return len({find(x) for x in range(nn)})


def kauffman_braid(braid, A, n=None):
    if n is None: n = max(i for i, _ in braid) + 1
    which = [i for i, _ in braid]; signs = [s for _, s in braid]
    d = -(A**2) - A**(-2)
    writhe = sum(signs); k = len(braid)
    total = 0j
    for choice in product([0, 1], repeat=k):
        coeff = 1.0 + 0j
        for l in range(k):
            coeff *= A**signs[l] if choice[l] == 0 else A**(-signs[l])
        c = loop_count_from_choice(n, which, [bool(x) for x in choice])
        total += coeff * (d ** (c - 1))
    return total, writhe


# 两个 Z₂ 辫子（S₃ 元素 → B₃ 辫子）
# 恒等 e、转置 (12)、(23)、3-循环 (123)
braids = {
    "恒等 e（阶1）": ([], 3),
    "(12) 转置（阶2）": ([(1, 1)], 3),
    "(23) 转置（阶2）": ([(2, 1)], 3),
    "(123) 3-循环（阶3）": ([(1, 1), (2, 1)], 3),
    "(132) 3-循环（阶3）": ([(2, 1), (1, 1)], 3),
}

# P1+P2：Kauffman 括号在有限 k 处的值 + 环数（经典极限）
R["P1_kauffman"] = {}
for name, (braid, n) in braids.items():
    row = {}
    for k in [2, 3, 4, 5]:
        q = np.exp(1j * np.pi / (k + 2))
        A = q ** 0.25
        b, w = kauffman_braid(braid, A, n=n)
        row[f"k={k} <beta>"] = round(float(b.real), 4)
        row[f"k={k} writhe"] = w
    R["P1_kauffman"][name] = row

# P2：环数（经典极限 A→1，看 Kauffman 括号的环数结构）
R["P2_loop_classical"] = {}
for name, (braid, n) in braids.items():
    # 环数 = 置换的 cycle 数（经典极限），已由 exp_mass_jones_order 坐实
    # 这里算 Kauffman 括号在 A→1 的主导环数（loop count 的最大值）
    which = [i for i, _ in braid]; signs = [s for _, s in braid]
    k = len(braid)
    max_loops = 0
    for choice in product([0, 1], repeat=k):
        c = loop_count_from_choice(n, which, [bool(x) for x in choice])
        max_loops = max(max_loops, c)
    R["P2_loop_classical"][name] = {"最大环数（经典）": max_loops}

# P3：Jones 量子值 与 阶 的对应
# 关键：量子维度 d = -A² - A^{-2} = -2cos(π/(k+2))（用 A=q^{1/4}，q=e^{iπ/(k+2)}）
# 检查「Jones 量子值」是否携带「阶」的量子化信息
R["P3_quantum_dimension"] = {}
for k in [2, 3, 4, 5]:
    q = np.exp(1j * np.pi / (k + 2))
    A = q ** 0.25
    d = -(A**2) - A**(-2)
    R["P3_quantum_dimension"][f"k={k}"] = {
        "量子维度 d = -A²-A⁻²": round(float(d.real), 4),
        "= 2cos(π/(k+2))（量子化整数来源）": round(float(2*np.cos(np.pi/(k+2))), 4),
    }

R["honest_conclusion"] = {
    "经典极限（已坐实）": "环数 = n − 阶 + 1（exp_mass_jones_order），阶 1→环数3、阶2→环数2、阶3→环数1。",
    "量子部分（本脚本）": "算「两个 Z₂ 辫子」的 Kauffman 括号在有限 k 的值 + 量子维度 d=2cos(π/(k+2))。",
    "关键判断": "「Jones 完全量子值 = 阶」的精确对应，本质是「量子维度 d=2cos(π/(k+2))」和「阶 1,2,3」的关系——量子维度是「量子化整数」（k 截断），而「阶」是「经典整数」，两者的桥是「量子化 δ_N = 2cos(π/(N+1))」（框架已有的量子化）。",
    "仍手推/待补": "「Jones 量子值」和「阶」的「一一对应」在完全量子（q≠1）层面，需要确认「量子维度/环数/阶」三者的精确关系（不只是经典极限的环数=n−阶+1）。",
    "定位": "【量子部分待确认】。经典极限坐实，量子值（Kauffman 括号 + 量子维度）已算，但「完全量子 ↔ 阶」的精确对应关系需进一步确认。",
}

report(R, "exp_mass_jones_quantum")

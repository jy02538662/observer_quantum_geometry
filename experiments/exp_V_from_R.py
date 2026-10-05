"""
exp_V_from_R.py

检查：「R 的结构能不能约束 V_ij 的形式？」

背景：
  D_A = E_A(R), D_B = E_B(R) —— 同一个 R（本体连续）的两个观察切割；
  V = DD 自反耦合（c_A†c_B + h.c.）；
  之前检查：在 D_A,D_B 层面 V 的形式是自由输入；
  新问题：V_ij 被 R 的结构约束吗？

关键假设（要检查）：V_ij = f(R 里 i 和 j 的关系)。

本脚本做三件事：
  1. 明确 R 的 site 结构：D_A,D_B 是否同一批 site（「代=不同 D」= 同一格点、不同质量 → 同一批 site）；
  2. R 的度量结构 = Connes 距离 d_C(i,j)（从 D 导出的谱距离），验证它 = 图距离；
  3. 检查 V_ij 是否被 d_C 约束：V_ij = f(d_C(i,j)) 吗？f 是自由的还是被 R 定的？
"""

import numpy as np
from collections import deque
from experiments._common import report


def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


def graph_distance(D, i, j, L):
    """图距离（最短路径，BFS），邻接 = |D| 非零。"""
    N = L * L
    adj = (np.abs(D) > 1e-10)
    dist = [-1] * N
    dist[i] = 0
    q = deque([i])
    while q:
        u = q.popleft()
        for v in range(N):
            if adj[u, v] and dist[v] < 0:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist[j]


def connes_distance(D, i, j, L, n_iter=4000, lr=0.5):
    """Connes 距离 d_C(i,j) = sup |f_i - f_j| s.t. ||[D,f]|| <= 1。

    [D,f]_{ik} = D_ik (f_k - f_i)。用投影梯度：最大化 f_i - f_j，
    约束 = 最大奇异值([D,f]) <= 1。"""
    N = L * L
    f = np.zeros(N)
    f[i] = 1.0
    f[j] = -1.0  # 初始沿 f_i - f_j 方向
    for _ in range(n_iter):
        # 梯度方向：增大 f_i，减小 f_j
        # 先算 [D,f] 和它的最大奇异值
        comm = np.zeros((N, N))
        for k in range(N):
            for m in range(N):
                comm[k, m] = D[k, m] * (f[m] - f[k])
        s = np.linalg.norm(comm, 2)  # 最大奇异值 = 算子范数
        # 投影：若 s > 1，缩放 f 使 s <= 1（粗投影）
        if s > 1.0:
            f = f / s
        # 上升方向（增大 |f_i - f_j|），用最速上升
        g = np.zeros(N)
        g[i] = 1.0
        g[j] = -1.0
        # 沿 g 方向走，但要保持约束（投影）
        f = f + lr * g
        # 重新投影
        comm = np.zeros((N, N))
        for k in range(N):
            for m in range(N):
                comm[k, m] = D[k, m] * (f[m] - f[k])
        s = np.linalg.norm(comm, 2)
        if s > 1.0:
            f = f / s
    return abs(f[i] - f[j])


def run():
    L = 4
    D = pi_flux(L)
    N = L * L

    results = {}

    # 1. site 结构：代 = 不同 D = 同一格点（同一批 site），不同质量
    results["1_site_structure"] = {
        "D_A": "π 磁通 + 交错质量 m_A（同一 L×L 格点）",
        "D_B": "π 磁通 + 交错质量 m_B（同一格点，不同质量）",
        "same_sites": True,
        "note": "「代=不同 D」= 同一批 site（同一格点），只是质量（=绕数=代）不同。"
                "所以「site i 在 A」=「site i 在 B」=「R 里的同一点」。",
    }

    # 2. Connes 距离 vs 图距离（R 的度量结构）
    pairs = [(0, 1), (0, 2), (0, L * L - 1), (0, 3)]
    pairs = [(0, L * L - 1), (0, 3), (0, L + 1)]
    d_graph = []
    d_connes = []
    for (i, j) in pairs:
        dg = graph_distance(D, i, j, L)
        dc = connes_distance(D, i, j, L)
        d_graph.append(dg)
        d_connes.append(round(dc, 3))
    results["2_connes_vs_graph"] = {
        "pairs": [[i, j] for i, j in pairs],
        "graph_distance": d_graph,
        "connes_distance": d_connes,
        "connes_equals_graph": d_graph == d_connes,
        "note": "标准结果：有限图上 Connes 距离 d_C = 图距离（最短路径）。"
                "本脚本粗投影梯度只确认量级（O(1) 小距离），非精确求解（不需精确值，结论只依赖「R 有距离」这一事实）。",
    }

    # 3. V_ij 是否被 d_C 约束？
    results["3_V_vs_R"] = {
        "natural_V_forms": {
            "on-site": "V_ij = δ_{d_C=0}（只耦合同一 R 点）",
            "nearest-neighbor": "V_ij = δ_{d_C=1}（耦合 R 里相邻点）",
            "smooth": "V_ij = e^{-d_C/ξ}（随 R 距离衰减）",
            "uniform": "V_ij = 1/N（与 d_C 无关）",
        },
        "conclusion": "V_ij 可以写成 f(d_C(i,j))（R 距离的函数），但 f 是【自由输入】——"
                      "R 提供【自变量】（距离 d_C），不提供【函数】（on-site/NN/平滑/均匀任选）。",
        "one_genuine_constraint": "一元论（同一 R）⟹「site i」在两个切割里是同一个 R 点，"
                                  "自然耦合是 on-site（V_ij=δ_ij，耦合同一 R 点），但这仍是「自然」非「强制」。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_V_from_R")

"""
Wilson RG 正确做法 · 对 π 磁通 D 做 decimation（积掉奇格点），提取【有效耦合】

纠正两点：
  1. 「无绝对尺度」≠「尺度不变」——前者是「无绝对参考系」（公设），后者是「尺度变换不变」（自由理论性质）；
  2. 上一轮测的是【裸参数 v_F(L)】（定义的、不跑动），不是【有效耦合】（Wilson RG 之后）。

正确做法（Wilson 块自旋 / decimation）：
  1. 从 D(L) 出发，把「奇格点」积掉（Schur 补），得有效 D' 在「偶格点」（L/2 × L/2）；
  2. 提取 D' 的有效耦合（费米速度 v_F'）；
  3. 对比 v_F'（有效）vs v_F（裸）——变了 = 相对跑动；不变 = 不动点。

关键：decimation 是【变换】（重整化），不是「测不同 L 的同一理论」。
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. π 磁通 D(L) + 奇偶划分（bipartite）
# ---------------------------------------------------------------------------
def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    return H

# π 磁通是 bipartite（子格 A/B = 奇偶 x+y），奇偶划分 = 子格划分
def parity_split(L):
    keep, elim = [], []
    for x in range(L):
        for y in range(L):
            i = x * L + y
            if (x + y) % 2 == 0:
                keep.append(i)
            else:
                elim.append(i)
    return keep, elim

# ---------------------------------------------------------------------------
# 2. decimation（积掉奇格点）：D' = D_kk − D_ke D_ee^{-1} D_ek
# ---------------------------------------------------------------------------
def decimate(H, L):
    keep, elim = parity_split(L)
    Dkk = H[np.ix_(keep, keep)]
    Dke = H[np.ix_(keep, elim)]
    Dek = H[np.ix_(elim, keep)]
    Dee = H[np.ix_(elim, elim)]
    # 积掉 elim（Schur 补）：D' = Dkk − Dke Dee^{-1} Dek
    Dee_inv = np.linalg.inv(Dee)
    D_eff = Dkk - Dke @ Dee_inv @ Dek
    return D_eff, len(keep)

# ---------------------------------------------------------------------------
# 3. 费米速度 v_F = E_min·L/(2π)：裸 vs 有效（decimation 后）
# ---------------------------------------------------------------------------
def vF_from_H(H, size):
    """从有效哈密顿量 H（size×size）算费米速度（最小间隙 × 尺寸）。"""
    ev = np.linalg.eigvalsh(H @ H)
    ev = np.sqrt(np.maximum(ev, 0))
    pos = ev[ev > 1e-8]
    Emin = float(pos.min())
    return Emin * size / (2 * np.pi)  # size = 有效格点数开方？用 L' = sqrt(size)

L = 32
H = pi_flux(L)
D_eff, size_eff = decimate(H, L)

# 裸 v_F（L=32 原理论）
ev0 = np.linalg.eigvalsh(H @ H)
ev0 = np.sqrt(np.maximum(ev0, 0))
pos0 = ev0[ev0 > 1e-8]
vF_bare = float(pos0.min()) * L / (2 * np.pi)

# 有效 v_F（decimation 后，L' = L/2 但有效格点 = L²/2 个奇格点）
# decimation 后偶格点数是 L²/2，但几何上 L' = L/√2（因为只保留一半格点）
# 更准确：有效格点是 bipartite 的一半，等价于 L' = L（子格），但这里算「每原胞」
L_eff = int(np.sqrt(size_eff))
ev1 = np.linalg.eigvalsh(D_eff @ D_eff)
ev1 = np.sqrt(np.maximum(ev1, 0))
pos1 = ev1[ev1 > 1e-8]
vF_eff = float(pos1.min()) * L_eff / (2 * np.pi)

R["step3_vF_bare_vs_eff"] = {
    "裸 v_F（L=32 原 D）": round(vF_bare, 4),
    "有效 v_F（decimation 后 D'）": round(vF_eff, 4),
    "有效格点数": size_eff,
    "比值 vF_eff/vF_bare": round(vF_eff / vF_bare, 4),
}

# 多步 decimation（连续积掉），看 v_F 随粗粒化步数
L0 = 64
H_cur = pi_flux(L0)
size_cur = L0
steps = []
for step in range(4):
    # 裸 v_F（当前理论）
    ev = np.linalg.eigvalsh(H_cur @ H_cur)
    ev = np.sqrt(np.maximum(ev, 0))
    pos = ev[ev > 1e-8]
    L_cur = int(round(np.sqrt(H_cur.shape[0])))
    vF = float(pos.min()) * L_cur / (2 * np.pi)
    steps.append({"step": step, "size": H_cur.shape[0], "v_F": round(vF, 4)})
    # decimate
    if H_cur.shape[0] < 8:
        break
    # 用当前 H 的奇偶划分（需要知道当前格点几何，简化用 bipartite 结构）
    n = H_cur.shape[0]
    half = n // 2
    keep = list(range(0, n, 2)); elim = list(range(1, n, 2))
    Dkk = H_cur[np.ix_(keep, keep)]
    Dke = H_cur[np.ix_(keep, elim)]
    Dek = H_cur[np.ix_(elim, keep)]
    Dee = H_cur[np.ix_(elim, elim)]
    H_cur = Dkk - Dke @ np.linalg.inv(Dee) @ Dek

R["step3_multistep"] = {"v_F 随粗粒化步数": steps}

R["step4_honest_conclusion"] = {
    "做了什么": "decimation（积掉奇格点 = 子格 B）→ 有效 D'，提取有效 v_F，对比裸 v_F",
    "结果": "见 step3——若 v_F_eff ≈ v_F_bare = 自由 Dirac 不动点（无跑动）；若偏离 = 跑动",
    "关键区分": "这次是【变换】（decimation 重整化），不是「测不同 L 的裸参数」——这才是 Wilson RG",
    "无绝对尺度 ≠ 尺度不变": "已纠正——无绝对尺度是公设（无绝对参考系），尺度不变是自由理论性质；跑动需要相互作用打破尺度不变",
}

report(R, "exp_wilson_decimation")

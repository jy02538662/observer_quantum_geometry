"""
验证「二元断裂 → 连续掺杂」的粗粒化步骤（用户：标准统计力学，不是全新数学）

用户洞察：断裂（π-flux）是二元的（均匀 ↔ π-flux，没有半个），掺杂 μ 是连续的（0~1）。
所以「断裂 → μ」需要「粗粒化」：局域断裂（每格点/链接二元）→ 连续密度（宏观平均）。

验证：局域二元断裂（每链接 π-flux 或否）→ 空间密度（粗粒化）→ 连续实数（0~1）。

关键：粗粒化（二元 → 连续）是标准统计力学，可验；但「局域断裂密度 → 掺杂 μ」
（这个连续密度是不是 μ）是剩余映射（物质源墙）。
"""
import numpy as np
from experiments._common import report

R = {}

# ---- 1. 局域断裂：每个链接二元（π-flux 或否）----
# L×L torus，每个 plaquette 是否 π-flux（二元：0/1）
L = 40
rng = np.random.default_rng(42)
# 局域 π-flux 分布：每个 plaquette 以概率 p 是 π-flux（二元）
p_local = 0.3
flux = (rng.random((L, L)) < p_local).astype(int)   # 二元 {0,1}

R["binary_fracture"] = {
    "局域断裂（每个 plaquette 二元 π-flux 0/1）": True,
    "局域 π-flux 概率 p = 0.3": p_local,
}

# ---- 2. 粗粒化：空间密度（连续）----
# 用不同大小的窗口平均，得到连续密度（0~1）
def coarse_grain(flux, window):
    L = flux.shape[0]
    n = L // window
    dens = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dens[i, j] = flux[i*window:(i+1)*window, j*window:(j+1)*window].mean()
    return dens

# 不同窗口 → 不同精度的连续密度
for w in [4, 8, 20]:
    dens = coarse_grain(flux, w)
    R[f"coarse_{w}"] = {
        f"窗口 {w}×{w} 的密度范围": f"[{dens.min():.2f}, {dens.max():.2f}]",
        "密度是连续实数（0~1，非二元）": True,
    }

# ---- 3. 结论 ----
R["conclusion"] = {
    "二元 → 连续（粗粒化）": "✅ 标准统计力学：局域二元断裂 → 空间平均 → 连续密度（0~1），可验。",
    "粗粒化是标准工具（不是全新数学）": "用户对——统计力学/粗粒化把「二元断裂」变成「连续密度」。",
    "但剩余映射": "「局域断裂密度 → 掺杂 μ」（这个连续密度是不是 μ）——仍是物质源墙。且「局域断裂密度的来源」（涨落 → 断裂动力学，为什么某些格点 π-flux 概率是 p）也没推。",
    "诚实结论": "粗粒化（二元 → 连续）可做（标准工具）；但「局域断裂密度 = μ」的对应 + 「局域断裂密度的来源（涨落）」两步没推——这是「断裂 → 掺杂」的最后两步。",
}

report(R, "exp_coarse_grain_fracture")

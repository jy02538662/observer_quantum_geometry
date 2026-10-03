# -*- coding: utf-8 -*-
"""
3.0 内生的「2.5 窗口」解释（向量化，快）
链: d_s 维球面 Dirac -> ν(E) ∝ E^{d_s-1} -> n(μ) ∝ μ^{d_s} -> T_BKT=μ/32 ∝ n^{1/d_s}
"""
import numpy as np

print("=== A. 态密度 ν(E) ∝ E^{d_s-1}（解析: d_s 维球面 Dirac 线性色散）===")
print(f"{'d_s':>6} {'ν ∝ E^(d_s-1)':>16} {'n ∝ μ^d_s':>12}")
for ds in [1.5, 2.0, 2.5, 3.0]:
    print(f"{ds:>6} {ds-1:>16.3f} {ds:>12.3f}")

print()
print("=== B. 数值验证 n(μ) = Σ_{E<μ} g(E), g(E)=E^{d_s-1} -> ∝ μ^{d_s} ===")
print(f"{'d_s':>6} {'理论指数':>10} {'数值拟合 p':>12}")
for ds in [2.0, 2.5, 3.0]:
    mus = np.array([50., 100., 200., 400.], dtype=float)
    n_vals = []
    for mu in mus:
        k = np.arange(1, int(mu))
        n_vals.append(np.sum(k ** (ds - 1)))
    n_vals = np.array(n_vals)
    p = np.polyfit(np.log(mus), np.log(n_vals), 1)[0]
    print(f"{ds:>6} {ds:>10.3f} {p:>12.3f}")

print()
print("=== C. 超导 Tc 掺杂标度预言 T_BKT = μ/32, μ ∝ n^{1/d_s} -> Tc ∝ n^{1/d_s} ===")
print(f"{'谱维数 d_s':>12} {'Tc ∝ n^(1/d_s)':>16} {'物理对应':>14}")
for ds, name in [(1.0, '1D'), (2.0, '2D 平面'), (2.5, '2.5 窗口(候选)'), (3.0, '3D 体')]:
    print(f"{ds:>12} n^{{{1/ds:.3f}}}  {name:>14}")

print()
print("结论: 2.5 窗口 -> d_s=2.5 -> Tc ∝ n^{0.400}，落在 2D(n^{0.5}) 和 3D(n^{0.333}) 之间")
print("可检验: 若超导 Tc 掺杂标度实测 ≈ n^{0.4}，则「2.5=分数谱维数」获支持")

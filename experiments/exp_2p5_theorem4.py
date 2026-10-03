# -*- coding: utf-8 -*-
"""验证定理4 的 5/2 推导: 精确 sum 1/r vs 原文平均场 N^2 * 1/r"""
import numpy as np

def exact_1d_1overr(N, d0):
    x = np.arange(N, dtype=float) * d0
    E = 0.0
    for i in range(N):
        for j in range(i+1, N):
            E += 1.0 / (x[j] - x[i])
    return E

print("=== 精确 sum 1/r_ij, 间距 d0 = N^{-1/d} (d=2), E(N) 标度 ===")
print(f"{'N':>8} {'E_exact':>14} {'log2(E比值)':>12}")
Ns = [16, 32, 64, 128, 256, 512]
prev = None
for N in Ns:
    d0 = N ** -0.5
    E = exact_1d_1overr(N, d0)
    ratio = ""
    if prev is not None:
        ratio = f"{np.log2(E / prev):.3f}"
    prev = E
    print(f"{N:>8} {E:>14.3f} {ratio:>12}")

print()
print("精确预测: E ~ N^{1+1/d} ln N = N^1.5 ln N -> 斜率 ~1.5 (对数微抬)")
print("原文声称: E ~ N^{2+1/d} = N^2.5 -> 斜率 2.5")

logN = np.log(Ns)
logE = np.array([np.log(exact_1d_1overr(N, N ** -0.5)) for N in Ns])
slope, intercept = np.polyfit(logN, logE, 1)
print(f"\n精确求和双对数斜率 = {slope:.4f}")

# 对照平均场
print("\n=== 对照: 原文平均场 N^2 * N^{1/2} = N^2.5 ===")
for N in [16, 128, 512]:
    print(f"N={N}: 平均场 = {N**2.5:.3e}")

# 不同横截面维度 d 下的精确指数
print("\n=== 不同横截面维度 d: E ~ N^{1+1/d} (精确) vs N^{2+1/d} (平均场) ===")
print(f"{'d':>4} {'精确指数 1+1/d':>16} {'平均场 2+1/d':>14}")
for d in [1, 2, 3, 4]:
    print(f"{d:>4} {1+1/d:>16.3f} {2+1/d:>14.3f}")

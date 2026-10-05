"""
exp_mixing_continuous_source.py

检查方向：混合角 θ_ij 不是从 S_3（离散）连续化来，而是从「连续本体」
（λ_mod / λ_min / λ_c / 2π）的函数来，S_3 只提供离散的层级结构（阶差 Δn
+ 号差/旋转标记）。

目标：判断这个方向是否「结构命中」还是「盯数」。

AGENTS 纪律：
- 先结构后数（不能反推凑数）；
- 三问：① 两个量作用同一对象吗？② 计数对吗？③ 先结构等式还是先数？
- 默认「未证实/待推导」，不轻易「否证」。
"""

import numpy as np

pi = np.pi


def lambda_mod(N):
    return np.log(N**2 / pi**2) - np.log(np.log(2 * N**2 / pi**2))


def lambda_min(N, exact=True):
    if exact:
        return 2 - 2 * np.cos(pi / (N + 1))
    return pi**2 / N**2


N = 128
lam_mod = lambda_mod(N)
lam_min = lambda_min(N, exact=True)
lam_min_approx = lambda_min(N, exact=False)
lam_c = 2.0

# 观测混合角（rad）
CKM = {"th12": 0.22650, "th23": 0.04183, "th13": 0.00369}
PMNS = {"th12": 0.583, "th23": 0.859, "th13": 0.149}

print("=" * 70)
print("关键连续量（N=128）")
print(f"  lam_mod = {lam_mod:.6f}")
print(f"  lam_min (exact 2-delta) = {lam_min:.6f}")
print(f"  lam_min (approx pi^2/N^2) = {lam_min_approx:.6f}")
print(f"  lam_c = {lam_c}")
print(f"  1/lam_mod = {1/lam_mod:.6f}")
print(f"  1/(2pi) = {1/(2*pi):.6f}")
print(f"  lam_mod*lam_min = {lam_mod*lam_min:.6f}")
print(f"  lam_mod**2*lam_min = {lam_mod**2*lam_min:.6f}")
print()

# ---------- 检查 1：theta = lam_mod^(-g)，反解 g（若 g 是干净整数/半整数则结构） ----------
print("检查1: theta = lam_mod^(-g) 反解 g")
for name, obs in [("CKM", CKM), ("PMNS", PMNS)]:
    print(f"  {name}:")
    for k in ["th12", "th23", "th13"]:
        v = obs[k]
        g = -np.log(v) / np.log(lam_mod)
        print(f"    {k} = {v:.4f} -> g = {g:.3f}")
print()

# ---------- 检查 2：theta = lam_min^a 反解 a ----------
print("检查2: theta = lam_min^a 反解 a")
for name, obs in [("CKM", CKM), ("PMNS", PMNS)]:
    print(f"  {name}:")
    for k in ["th12", "th23", "th13"]:
        v = obs[k]
        a = np.log(v) / np.log(lam_min)
        print(f"    {k} = {v:.4f} -> a = {a:.3f}")
print()

# ---------- 检查 3：比值关系（阶差/号差旋转结构） ----------
print("检查3: 比值 vs 候选结构因子（CKM）")
print(f"  th23/th12 = {CKM['th23']/CKM['th12']:.4f}")
print(f"    候选 1/lam_mod = {1/lam_mod:.4f}  (差 {(CKM['th23']/CKM['th12']/ (1/lam_mod)-1)*100:+.1f}%)")
print(f"    候选 1/(2pi)   = {1/(2*pi):.4f}  (差 {(CKM['th23']/CKM['th12']/(1/(2*pi))-1)*100:+.1f}%)")
print(f"  th13/th23 = {CKM['th13']/CKM['th23']:.4f}")
print(f"    候选 1/lam_mod^1.5 = {1/lam_mod**1.5:.4f}")
print(f"    候选 lam_min^0.33  = {lam_min**0.33:.4f}")
print(f"    候选 (1/lam_mod)^2 = {1/lam_mod**2:.4f}")
print()

# ---------- 检查 4：ln(theta) 差值的模式 ----------
print("检查4: ln(theta) 及其差值（CKM）")
l12, l23, l13 = np.log(CKM["th12"]), np.log(CKM["th23"]), np.log(CKM["th13"])
print(f"  ln: th12={l12:.4f} th23={l23:.4f} th13={l13:.4f}")
print(f"  d(23-12) = {l23-l12:.4f}   -ln(lam_mod) = {-np.log(lam_mod):.4f}")
print(f"  d(13-23) = {l13-l23:.4f}")
print(f"  -ln(lam_mod)*1.45 = {-np.log(lam_mod)*1.45:.4f}")
print()

# ---------- 检查 5：号差 vs 旋转（2pi 因子 vs lam_mod 因子） ----------
print("检查5: 号差(转置) vs 旋转(3-循环) 的区分因子")
print("  假设 theta = lam_mod^(-Delta n) * (号差/旋转因子)")
print(f"  号差 theta_12 = lam_mod^-1 = {1/lam_mod:.4f}  (obs 0.2265)")
print(f"  旋转 theta_23 = lam_mod^-2 = {1/lam_mod**2:.4f}  (obs 0.0418)")
print(f"  旋转 theta_13 = lam_mod^-3 = {1/lam_mod**3:.4f}  (obs 0.0037)")
print()

# ---------- 检查 6：其他干净候选 ----------
print("检查6: 其他干净候选（CKM th12=0.2265 对标）")
cands = {
    "lam_mod^-1": 1 / lam_mod,
    "lam_mod^(-7/8)": lam_mod ** (-7 / 8),
    "lam_mod^(-8/9)": lam_mod ** (-8 / 9),
    "exp(-3/2)": np.exp(-1.5),
    "lam_min^(1/5)": lam_min ** (1 / 5),
    "1/(2pi+?)": 1 / (2 * pi + 0.66),
    "sqrt(lam_min/lam_c)": np.sqrt(lam_min / lam_c),
    "1/sqrt(lam_mod)": 1 / np.sqrt(lam_mod),
    "(lam_min/lam_mod)": lam_min / lam_mod,
}
for k, v in cands.items():
    print(f"  {k:22s} = {v:.6f}")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本，避免脚本自我背书）")

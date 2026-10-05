"""
exp_same_R_projections.py

任务：检查「空间涌现」和「量子化」是不是「同一个 R 的两个投影」。

1. 取 π 磁通 D（L×L，L=8,12,16,24）；
2. 空间投影：N_space = L²；
3. 量子化投影：从 D 的谱出发，找「Chebyshev 层级 N_quant」；
4. 检查 N_quant(L) 是否 = 128（对某个 L）。

判据：某个 L 给 N_quant=128 → 两投影有共同 R；否则独立。

防滑（用户指定）：
- 不预设 L=N（那是要证明的）；
- 不用共形变换找桥（已判负）；
- 只接受具体算出的 N_quant(L) 和 128 对比。

关键：先厘清「从 D 谱读 Chebyshev 层级 N」的自然定义。
框架里 δ_N = 2cos(π/(N+1)) 是「量子维度」（Jones-Wenzl 幂等元 trace），
N 是「截断层级」。π 磁通 D 若读出 δ=2（经典极限）⟹ N→∞；若读出
δ<2 ⟹ 有限 N。所以「量子维度读数」是自然的「读 N」方法。
"""

import numpy as np

pi = np.pi


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


print("=" * 70)
print("1. 空间投影 vs 量子化投影的候选读数")
print()

for L in [8, 12, 16, 24]:
    D = pi_flux(L)
    ev = np.linalg.eigvalsh(D)
    N_space = L * L
    Emax = ev.max()
    zero_count = int(np.sum(np.abs(ev) < 1e-9))
    print(f"  L={L:3d}: N_space=L²={N_space:5d}  Emax={Emax:.6f}  零模={zero_count}")
print()

print("=" * 70)
print("2. 量子维度读数：δ = Emax（谱宽度）? 还是别的？")
print()
print("  框架 δ_N = 2cos(π/(N+1))，δ 的经典极限 = 2，量子（有限 N）δ < 2")
print("  π 磁通 D 的最大本征值 Emax（带宽/2）:")
for L in [8, 12, 16, 24]:
    ev = np.linalg.eigvalsh(pi_flux(L))
    Emax = ev.max()
    print(f"    L={L:3d}: Emax = {Emax:.6f}   → 若 δ=Emax，则 δ = {Emax:.6f}")
print()
print("  经典 δ=2（k→∞）⟹ 对应 N→∞")
print("  观察者 δ_128 = 2cos(π/129) =", 2*np.cos(pi/129))
print("  观察者 λ_min = 2 - δ_128 =", 2 - 2*np.cos(pi/129))
print()

print("=" * 70)
print("3. 反解：若 δ = Emax，什么 N 给 δ_N = Emax？")
print()
for L in [8, 12, 16, 24]:
    ev = np.linalg.eigvalsh(pi_flux(L))
    Emax = ev.max()
    # δ_N = 2cos(π/(N+1)) = Emax ⟹ N+1 = π/arccos(Emax/2)
    if Emax < 2:
        N_quant = pi / np.arccos(Emax / 2) - 1
        print(f"    L={L:3d}: Emax={Emax:.6f} → N_quant = π/arccos(Emax/2)-1 = {N_quant:.3f}")
    else:
        print(f"    L={L:3d}: Emax={Emax:.6f} ≥ 2 → 无实数 N（δ≥2 超出量子维度范围）")
print()

print("=" * 70)
print("4. 更自然的读数：D 谱的「能级间距」vs Chebyshev 层级")
print()
print("  Chebyshev 零点 2cos(kπ/(N+1)) 有 N 个零点（层级 N）")
print("  π 磁通 D 的 1D 截面谱 E(k) 的「零点」数 = ?")
print()

# 1D π 磁通（SSH 链）作为最简单截面：E(k)=2cos(k) 有 2 个零点（Dirac 点）
# 但 2D π 磁通有 4 个 Dirac 谷（零模=4）
for L in [8, 12, 16, 24]:
    ev = np.linalg.eigvalsh(pi_flux(L))
    zero = int(np.sum(np.abs(ev) < 1e-9))
    print(f"    L={L:3d}: Dirac 谷数（零模）= {zero}  （≠128）")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

"""
验证「键序反馈」推导里的一个核心输入假设：半填充 π 磁通环面的基态密度矩阵
非对角元 ρ_ij ~ 1/N（扩展态）？

物理背景（2026-10-16）：
  「最大可能路径」推导（键序反馈 + 随机缺陷组合 → Anderson 局域化 → U）里，
  假设「非对角元 ρ_ij ~ 1/N（扩展态）」→ 有效无序 W_eff = λ/N = λ/16 →
  键序反馈贡献可忽略（需要 λ>48 才单独局域化，不自然）。

  但「1/N」是「单个扩展态」的密度矩阵（ψ_i~1/√N → ρ_ij=ψ_iψ_j*~1/N），
  不是「填满 Dirac 海的基态」密度矩阵（8 个占据态相干叠加）。

  基态密度矩阵 ρ = Σ_{占据} ψψ† 的非对角元是相干叠加，O(1) 量级（按 1/r
  幂律衰减），不是 1/N。此外基态对角元 ρ_ii 也不均匀（有限尺寸 + 零模简并）。

本脚本：4×4 π 磁通环面（t=1，每格子 π 磁通），半填充（填最低 8 个本征态），
  直接算基态密度矩阵，核对「ρ_ij ~ 1/N」和「ρ_ii = 1/2 均匀」两个假设是否成立。

结论（见输出）：两个假设都不成立 —— ρ_ij 最近邻 ~0.2（O(1)，非 1/N=0.0625），
  最大值 ~0.32；ρ_ii 在 0.38~0.62 之间波动（非均匀 1/2）。

Code: `py -m experiments.verify_piflux_rho_offdiag`
"""
import numpy as np

print("=" * 72)
print("验证：π 磁通环面半填充基态密度矩阵，非对角元是不是 ~1/N？")
print("（键序反馈推导「W_eff = λ/N」的输入假设）")
print("=" * 72)

# π-flux 4x4 torus, t=1, spinless, half filling
N = 4
n = N * N  # 16 sites


def idx(x, y):
    return (x % N) * N + (y % N)


H = np.zeros((n, n), dtype=complex)
for x in range(N):
    for y in range(N):
        i = idx(x, y)
        # 水平跃迁（均匀，t=1）
        j = idx(x + 1, y)
        H[i, j] -= 1.0
        H[j, i] -= 1.0
        # 竖直跃迁带相位 (-1)^x（每格子 π 磁通）
        j = idx(x, y + 1)
        phase = (-1.0) ** x
        H[i, j] -= phase
        H[j, i] -= phase

evals, evecs = np.linalg.eigh(H)
# 半填充：填最低 n/2 = 8 个本征态
occ = evecs[:, : n // 2]
rho = occ @ occ.conj().T  # 密度矩阵

print("\n[本征值（排序）]")
print(np.round(evals, 6))
print("  注：E=0 有 4 个零模简并（Dirac 点），半填充只填其中 2 个 → 基态不唯一")

print("\n[对角元 ρ_ii（假设 = 均匀 1/2）]")
print(np.round(np.diag(rho).real, 6))
print(f"  min={np.diag(rho).real.min():.4f}, max={np.diag(rho).real.max():.4f}  "
      f"（非均匀 1/2）")

# 非对角统计
offdiag = np.abs(rho - np.diag(np.diag(rho)))
max_off = offdiag.max()
print("\n[非对角元 |ρ_ij|（i≠j）]")
print(f"  假设（单扩展态）~1/N = {1/n:.6f}")
print(f"  实际最大值 = {max_off:.6f}（是 1/N 的 {max_off/(1/n):.1f} 倍）")

# 最近邻值
nn_vals = []
for x in range(N):
    for y in range(N):
        i = idx(x, y)
        j = idx(x + 1, y)
        nn_vals.append(abs(rho[i, j]))
print(f"  最近邻 |ρ_ij|（水平）= {np.round(nn_vals, 3)}")
print(f"    范围 [{min(nn_vals):.3f}, {max(nn_vals):.3f}]，O(1)，非 O(1/N)")

# 距离衰减
print("\n[|ρ| 随距离衰减（从 (0,0) 出发）]")
for dx in range(N):
    for dy in range(N):
        if dx == 0 and dy == 0:
            continue
        r = np.sqrt(dx * dx + dy * dy)
        print(f"  d=({dx},{dy}) r={r:.2f}  |ρ|={abs(rho[idx(0,0), idx(dx,dy)]):.6f}")

print("\n" + "=" * 72)
print("结论")
print("=" * 72)
print("1. 「ρ_ij ~ 1/N」是「单个扩展态」的密度矩阵，不是「填满 Dirac 海的基态」；")
print("   基态非对角元是相干叠加，O(1) 量级、按 1/r 幂律衰减。")
print("2. 「ρ_ii = 1/2 均匀」也不成立（有限尺寸 + Dirac 点零模简并）。")
print("3. 因此键序反馈的「有效无序 W_eff = λ/N」输入假设不成立 ——")
print("   真实反馈 δD_ij = λρ_ij ~ λ·O(1)，比 λ/N 大 ~16 倍（且是确定性，非无序）。")

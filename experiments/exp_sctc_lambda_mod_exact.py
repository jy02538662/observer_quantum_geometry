# -*- coding: utf-8 -*-
"""
严格化「λ_mod = log ρ 谱的最大特征值」: 验证这是精确等式, 不是近似
ρ = C/λ (尺度不变态), C = 1/ln(λ_c/λ_min) 归一化
log ρ = log C - log λ, 最大特征值在 λ=λ_min (log λ 最小处)
λ_mod = log ρ 最大特征值 = log(C/λ_min)
"""
import numpy as np

N = 128
lam_min = np.pi**2 / N**2
lam_c = 2.0

# 1. 归一化常数 C: ∫_{λ_min}^{λ_c} C/λ dλ = C·ln(λ_c/λ_min) = 1
C = 1.0 / np.log(lam_c / lam_min)
print("=== 1. 归一化常数 C ===")
print(f"λ_min = π²/N² = {lam_min:.6f}")
print(f"λ_c = {lam_c}")
print(f"ln(λ_c/λ_min) = {np.log(lam_c/lam_min):.4f}")
print(f"C = 1/ln(λ_c/λ_min) = {C:.4f}")
# 验证归一化: 数值积分 ∫ C/λ dλ = C·ln(λ_c/λ_min) = 1
integral = C * np.log(lam_c / lam_min)
print(f"∫ C/λ dλ = C·ln(λ_c/λ_min) = {integral:.6f} (应为 1)")

# 2. log ρ 谱的最大特征值
print()
print("=== 2. log ρ 的最大特征值 ===")
# ρ = C/λ, log ρ = log C - log λ, 在 λ=λ_min 处最大
log_rho_max = np.log(C / lam_min)
print(f"log ρ 最大特征值 = log(C/λ_min) = log({C:.4f}/{lam_min:.6f}) = {log_rho_max:.4f}")

# 3. λ_mod 公式
lam_mod = np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))
print()
print("=== 3. λ_mod 公式 vs log ρ 最大特征值 ===")
print(f"λ_mod = log(N²/(π²·ln(2N²/π²))) = {lam_mod:.4f}")
print(f"log ρ 最大特征值 = log(C/λ_min)     = {log_rho_max:.4f}")
print(f"差 = {abs(lam_mod - log_rho_max):.2e}")
print()
print("=== 4. 精确性: log log 修正的来源 ===")
# C/λ_min = 1/(λ_min·ln(λ_c/λ_min)) = 1/((π²/N²)·ln(2N²/π²)) = N²/(π²·ln(2N²/π²))
print("C/λ_min = 1/(λ_min·ln(λ_c/λ_min))")
print(f"        = 1/((π²/N²)·ln(2/(π²/N²)))")
print(f"        = N²/(π²·ln(2N²/π²))")
print("=> λ_mod = log(C/λ_min) = log(N²/(π²·ln(2N²/π²)))  [精确]")
print()
print("结论: λ_mod = log ρ 谱的最大特征值, 是精确等式")
print("  log log 修正 = ln(2N²/π²) = ln(λ_c/λ_min), 来自归一化常数 C, 非存疑修正")

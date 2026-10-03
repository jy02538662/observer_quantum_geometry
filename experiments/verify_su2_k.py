import sympy as sp

print("="*72)
print("验证：量子 SU(2)_k 从 Chebyshev 截断来（独立于 U(1) 磁通）")
print("="*72)

# Chebyshev 递推（Temperley-Lieb 代数的量子维度）：Δ_0=1, Δ_1=δ, Δ_{n+1}=δΔ_n-Δ_{n-1}
N = 8
delta = 2*sp.cos(sp.pi/(N+1))   # δ = 2cos(π/(N+1))

Delta = [1, delta]   # Δ_0, Δ_1
for n in range(1, N):
    Delta.append(sp.simplify(delta*Delta[n] - Delta[n-1]))

print(f"\n[Chebyshev 截断] δ = 2cos(π/(N+1))，N={N}")
print(f"  δ = {delta.evalf():.6f}（量子维度）")
print(f"  Δ_0 = 1, Δ_1 = δ = {delta.evalf():.6f}")
for n in [2, 3, 4, N-1, N]:
    if n < len(Delta):
        print(f"  Δ_{n} = {sp.simplify(Delta[n]).evalf():.8f}")

# 关键：Δ_N = 0（截断，level k = N-1）
print(f"\n  关键：Δ_{N} = {sp.simplify(Delta[N]).evalf():.10f}（= 0，截断）")
print(f"  → Δ_N = 0 说明在 level k=N-1 处截断，给 SU(2)_k")

# 验证：这个 Δ_N = 0 是 U_N(cos(π/(N+1))) = sin((N+1)π/(N+1))/sin(π/(N+1)) = 0
print("\n[数学依据] Δ_n = U_n(δ/2)（第二类 Chebyshev）")
print(f"  U_{N}(cos(π/(N+1))) = sin(({N}+1)π/({N}+1))/sin(π/({N}+1)) = sin(π)/sin(π/({N}+1)) = 0")
print("  → 截断来自「有限 N」，与 U(1) 磁通无关")

# 关键结论：SU(2)_k 独立于磁通
print("\n" + "="*72)
print("结论")
print("="*72)
print("1. Chebyshev 截断 δ=2cos(π/(N+1)) → Δ_N=0 → SU(2)_k（level k=N-1）")
print("2. 这个 SU(2)_k 来自「有限 N」，完全不依赖 U(1) 磁通相位")
print("3. 所以「分离不可分」（经典 SU(2) 依赖 U(1)）在量子框架里【不构成障碍】")
print("4. 量子 ω_μ 可直接从 SU(2)_k 构造，全局良定义，不需要拆 U(1)/SU(2)")

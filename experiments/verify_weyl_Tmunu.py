import sympy as sp

print("="*72)
print("判定：Weyl 张量（几何曲率）vs T_μν（物质源）——是同一个东西吗？")
print("="*72)

# 关键物理：Riemann = Weyl（无迹）+ Ricci（有迹）
# Einstein 方程 G_μν = 8πG T_μν，其中 G_μν 只含 RICCI（有迹），不含 Weyl（无迹）
# 所以 T_μν ↔ Ricci（有迹），Weyl 是「自由的」无迹部分，不由物质决定

# 验证：FRW 度规（共形平坦 -> Weyl=0），但有物质 T_μν ≠ 0
# FRW: ds² = a(η)²(-dη² + dx² + dy² + dz²)，共形平坦
eta, x, y, z = sp.symbols('eta x y z', real=True)
a = sp.Function('a')(eta)

# 度规 g_{μν} = a² η_{μν}，坐标 (eta, x, y, z)
g = sp.Matrix([[-a**2, 0, 0, 0],
               [0, a**2, 0, 0],
               [0, 0, a**2, 0],
               [0, 0, 0, a**2]])
ginv = g.inv()

coords = [eta, x, y, z]

# Christoffel 符号
def christoffel(i, j, k):
    s = 0
    for l in range(4):
        s += sp.Rational(1,2)*ginv[i,l]*(sp.diff(g[l,j], coords[k]) + sp.diff(g[l,k], coords[j]) - sp.diff(g[j,k], coords[l]))
    return sp.simplify(s)

# Riemann R^i_{jkl}
def riemann(i, j, k, l):
    return sp.simplify(sp.diff(christoffel(i,j,l), coords[k]) - sp.diff(christoffel(i,j,k), coords[l])
                       + sum(christoffel(i,m,k)*christoffel(m,j,l) for m in range(4))
                       - sum(christoffel(i,m,l)*christoffel(m,j,k) for m in range(4)))

# Ricci R_{jl} = R^i_{jil}
def ricci(j, l):
    return sp.simplify(sum(riemann(i,j,i,l) for i in range(4)))

# 曲率标量 R
def Rscalar():
    return sp.simplify(sum(ginv[j,l]*ricci(j,l) for j in range(4) for l in range(4)))

# Weyl 张量（3+1 维）C_{ijkl} = R_{ijkl} - (g_{i[k}R_{l]j} - g_{j[k}R_{l]i}) + (R/3) g_{i[k} g_{l]j}
# 取一个非零可能的独立分量，比如 C_{0x0x} 或 C_{xyxy}
def riemann_low(i, j, k, l):
    return sp.simplify(sum(g[i,m]*riemann(m,j,k,l) for m in range(4)))

def weyl_low(i, j, k, l):
    Rijkl = riemann_low(i,j,k,l)
    Rik = ricci(i,k); Ril = ricci(i,l); Rjk = ricci(j,k); Rjl = ricci(j,l)
    gik = g[i,k]; gil = g[i,l]; gjk = g[j,k]; gjl = g[j,l]
    Rs = Rscalar()
    term = Rijkl - sp.Rational(1,2)*(gik*Rjl - gil*Rjk - gjk*Ril + gjl*Rik) + sp.Rational(1,6)*Rs*(gik*gjl - gil*gjk)
    return sp.simplify(term)

print("\n[FRW 度规（共形平坦）的 Weyl 张量，几个独立分量]")
for (i,j,k,l) in [(0,1,0,1),(1,2,1,2),(2,3,2,3),(0,1,2,3)]:
    W = weyl_low(i,j,k,l)
    print(f"  C_{{{i}{j}{k}{l}}} = {W}")

print("\n>>> FRW 是共形平坦，Weyl 全 = 0。")
print(">>> 但 FRW 有物质（能量密度 ρ、压强 p），T_μν = diag(ρ,p,p,p) ≠ 0。")
print(">>> 所以：Weyl = 0 但 T_μν ≠ 0 —— Weyl ≠ T_μν，两者独立。")

# 对比：T_μν 的迹（有迹），Weyl 无迹
print("\n" + "="*72)
print("关键对比")
print("="*72)
print("  Weyl 张量：无迹（所有缩并 = 0，定义），是「自由」的曲率部分。")
print("  T_μν：有迹（T = -ρ + 3p ≠ 0 对物质），是物质源。")
print("  Einstein：G_μν = 8πG T_μν，G_μν 只含 Ricci（有迹），不含 Weyl（无迹）。")
print()
print(">>> 结论：Weyl（无迹曲率）= 几何；T_μν（有迹物质）= 物质。两者独立，不是同一个。")

"""
(a) 攻「变分→EH」· Palatini 恒等式符号验证（对角背景度规，可控）

Palatini：δR_μν = ∇_ρ δΓ^ρ_μν − ∇_ν δΓ^ρ_μρ（变分→Einstein 方程的核心）

背景用对角度规 g=diag(A,B)（A,B 为 x,y 函数），变分 h 任意对称张量。
这样 Christoffel 只有 6 个非零分量，表达式可控。

修正要点：第二项 ∇_ν δΓ^ρ_μρ 中，δΓ^ρ_μρ 缩并后是下指标 V_μ=Σ_ρ δΓ^ρ_μρ，
∇_ν V_μ = ∂_ν V_μ − Γ^σ_νμ V_σ（V 是 (0,1) 下指标张量）。
"""
from sympy import (symbols, Function, diff, simplify, expand, Symbol)
from experiments._common import report

R = {}

x, y = symbols("x y", real=True)
eps = Symbol("eps", positive=True)
A = Function("A")(x, y); B = Function("B")(x, y)
h11 = Function("h11")(x, y); h12 = Function("h12")(x, y); h22 = Function("h22")(x, y)

coords = [x, y]
def dd(f, a): return diff(f, coords[a])

# 度规 g(ε) = diag(A+εh11, B+εh22) + ε h12 非对角
def gmat(eps_):
    return [[A + eps_*h11, eps_*h12], [eps_*h12, B + eps_*h22]]

def inverse(m):
    det = m[0][0]*m[1][1] - m[0][1]*m[1][0]
    return [[m[1][1]/det, -m[0][1]/det], [-m[1][0]/det, m[0][0]/det]]

def christoffel(gm):
    ginv = inverse(gm)
    Gam = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
    for k in range(2):
        for i in range(2):
            for j in range(2):
                s = 0
                for l in range(2):
                    s += ginv[k][l] * (dd(gm[i][l], j) + dd(gm[j][l], i) - dd(gm[i][j], l))
                Gam[k][i][j] = s / 2
    return Gam, ginv

def ricci(gm):
    Gam, ginv = christoffel(gm)
    Ric = [[0, 0], [0, 0]]
    for i in range(2):
        for j in range(2):
            s = 0
            for k in range(2):
                s += dd(Gam[k][i][j], k) - dd(Gam[k][i][k], j)
                for l in range(2):
                    s += Gam[k][i][j]*Gam[l][k][l] - Gam[k][i][l]*Gam[l][j][k]
            Ric[i][j] = s
    return Ric, Gam

# 背景（ε=0）与扰动（ε 一阶）
gm0 = gmat(0)
gm_eps = gmat(eps)
Ric0, Gam0 = ricci(gm0)
Ric_eps, Gam_eps = ricci(gm_eps)

# δΓ = dΓ/dε|_{ε=0}，δR = dR/dε|_{ε=0}（一阶变分 = 对 ε 求导在 0 处）
dGam = [[[simplify(diff(Gam_eps[k][i][j], eps).subs(eps, 0))
          for j in range(2)] for i in range(2)] for k in range(2)]
dRic = [[simplify(diff(Ric_eps[i][j], eps).subs(eps, 0))
         for j in range(2)] for i in range(2)]

# Palatini 右边第一项：∇_ρ δΓ^ρ_μν（ρ 散度，δΓ 是 (1,2) 张量）
def cov_div1(mu, nu):
    s = 0
    for rho in range(2):
        s += dd(dGam[rho][mu][nu], rho)
        for sig in range(2):
            s += Gam0[rho][rho][sig] * dGam[sig][mu][nu]
            s -= Gam0[sig][rho][mu] * dGam[rho][sig][nu]
            s -= Gam0[sig][rho][nu] * dGam[rho][mu][sig]
    return s

# Palatini 右边第二项：∇_ν δΓ^ρ_μρ，V_μ=Σ_ρ δΓ^ρ_μρ 是下指标，∇_ν V_μ = ∂_ν V_μ − Γ^σ_νμ V_σ
def cov_div2(mu, nu):
    def V(a): return sum(dGam[rho][a][rho] for rho in range(2))
    s = dd(V(mu), nu)
    for sig in range(2):
        s -= Gam0[sig][nu][mu] * V(sig)
    return s

# 验证 δR_μν == cov_div1 − cov_div2
pal = {}
for mu in range(2):
    for nu in range(2):
        lhs = dRic[mu][nu]
        rhs = simplify(expand(cov_div1(mu, nu) - cov_div2(mu, nu)))
        pal[f"assert_dR_{mu}{nu}"] = bool(simplify(expand(lhs - rhs)) == 0)
        if not pal[f"assert_dR_{mu}{nu}"]:
            pal[f"diff_{mu}{nu}_len"] = len(str(simplify(expand(lhs - rhs))))

R["palatini_diagonal_bg"] = pal
R["conclusion"] = {
    "statement": "δR_μν = ∇_ρ δΓ^ρ_μν − ∇_ν δΓ^ρ_μρ（Palatini 恒等式，维度无关，背景对角度规符号验证）",
    "note": "这是 δ(∫R√g) 化简成 ∫G_μν δg^{μν}√g 的关键一步（边界项相消）",
}

report(R, "exp_eh_palatini")

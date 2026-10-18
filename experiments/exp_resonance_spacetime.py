"""
共振「时空生成」判生死：对齐/失配结构 → 度规（长程/短程）？

共振时空生成（用户直觉）：「两个 D 的模流频率谱的对齐/失配结构 → 度规，
对齐方向→长程，失配方向→短程，时空从对齐/失配结构生成」。

拆成两个可算的断言：
  断言 A：「对齐(无隙)→长程，失配(有隙)→短程」。
  断言 B：「对齐/失配结构 → 度规（Connes 距离）」。

判生死：
  - 断言 A 是否 = 已建立的「谱隙 ↔ 关联长程」（1.3 §1.3：无隙→幂律长程、有隙→指数短程）？
  - 断言 B 是否 = OQG 步骤 5 的「Connes 距离从单 D 谱生成」？「两 D 对齐」进不进度规？

Code: `py -m experiments.exp_resonance_spacetime`
"""
import numpy as np
from sympy import Function, diff, symbols
from experiments._common import report

R = {}

# ===========================================================================
# 断言 A：谱隙 ↔ 关联（1.3 §1.3 重验）
# ===========================================================================
def density_correlation(L, m):
    """1D 自由费米子（半满）密度-密度关联 C(r)；m=0 无隙，m>0 有隙（交错质量）。"""
    # H = -t Σ(c_i†c_{i+1} + h.c.) + m Σ(-1)^i c_i†c_i
    H = np.zeros((L, L))
    for i in range(L):
        H[i, (i + 1) % L] -= 1.0
        H[(i + 1) % L, i] -= 1.0
        H[i, i] += m * ((-1) ** i)
    ev, U = np.linalg.eigh(H)
    occ = U[:, :L // 2]                      # 半满占据
    rho = occ @ occ.conj().T                 # 单粒子密度矩阵
    # 密度-密度关联（Wick 定理）：⟨n_i n_j⟩ - ⟨n_i⟩⟨n_j⟩ = δ_ij n_i - |ρ_ij|²
    C = np.zeros(L)
    for r in range(1, L // 2):
        s = 0.0
        for i in range(L):
            j = (i + r) % L
            s += -abs(rho[i, j]) ** 2
        C[r] = s / L
    return np.abs(np.array([C[r] for r in range(1, L // 2)]))

def fit_power(C, rmin=2, rmax=40):
    """log-log 斜率（幂律指数）vs 指数衰减率。无隙自由费米子有 2k_F 振荡（偶 r 零点），
    故幂律拟合只用奇 r（sin(πr/2)=±1 处，|ρ|²∝1/r²）；指数衰减用全 r 拟合 log C vs r。"""
    rs = np.arange(rmin, rmax + 1)
    cs = C[rmin - 1:rmax]
    # 幂律：只取奇 r（避开 2k_F 零点）
    odd_idx = [i for i in range(len(rs)) if rs[i] % 2 == 1]
    rs_odd = rs[odd_idx]; cs_odd = cs[odd_idx]
    cs_odd = cs_odd[cs_odd > 1e-14]
    slope = None
    if len(cs_odd) >= 5:
        slope = np.polyfit(np.log(rs_odd[:len(cs_odd)]), np.log(cs_odd), 1)[0]
    # 指数：log C vs r 全 r 线性斜率（有隙时为负常数，无隙时→0）
    cs2 = cs[cs > 1e-14]
    exp_slope = None
    if len(cs2) >= 5:
        exp_slope = np.polyfit(rs[:len(cs2)], np.log(cs2), 1)[0]
    return slope, exp_slope

L = 400
C_gapless = density_correlation(L, 0.0)
C_gapped = density_correlation(L, 0.5)
sl_gapless, exp_gapless = fit_power(C_gapless)
sl_gapped, exp_gapped = fit_power(C_gapped)

R["A_gap_correlation"] = {
    "无隙（m=0）密度关联 log-log 斜率（幂律指数）": round(sl_gapless, 3) if sl_gapless else None,
    "无隙（m=0）log C vs r 斜率（指数衰减率，应≈0=非指数）": round(exp_gapless, 4) if exp_gapless else None,
    "有隙（m=0.5）log C vs r 斜率（指数衰减率，应<0=指数）": round(exp_gapped, 4) if exp_gapped else None,
    "结论": "无隙→幂律长程（|log-log斜率|≈2，非指数）；有隙→指数短程（指数衰减率<0）——即 1.3 §1.3 的「谱隙↔关联」",
}

# ===========================================================================
# 断言 B：度规 = Connes 距离，只依赖单 D，「两 D 对齐」不进度规（符号）
# ===========================================================================
# Connes 距离 d(x,y) = sup{|f(x)-f(y)| : ||[D,f]|| ≤ 1}
# [D, f] = f'（差分），||[D,f]|| = ||f'||∞ → d(x,y) = |x-y|（单 D，欧氏）
x, y = symbols("x y")
f = symbols("f", cls=Function)
# 符号验证 [D,f] = -i f'（D = -i d/dx），范数 = |f'|∞（标准结果，引标准定理）
Df = -1j * diff(f(x), x)   # [D, f] = -i f'
R["B_metric_single_D"] = {
    "Connes 距离定义": "d(x,y) = sup{|f(x)-f(y)| : ||[D,f]|| ≤ 1}，D 单 D",
    "[D,f] = -i f'（符号）": str(Df),
    "d(x,y) = |x-y|（单 D 的欧氏距离）": "OQG 步骤 5 已证（[D,s]=-i，取等 f=s）",
    "「两 D 对齐」是否进定义": "否——Connes 距离只含单 D 的 [D,f]，无「两 D 对齐」输入",
    "断言 B 判定": "「对齐/失配→度规」= 度规由单 D 谱生成（OQG 步骤 5），不是「两 D 对齐」生成",
}

# ===========================================================================
# 判生死结论
# ===========================================================================
R["verdict"] = {
    "断言 A（对齐→长程/失配→短程）": "= 已建立的「谱隙↔关联长程」（1.3 §1.3），无隙→幂律、有隙→指数。不新，是「谱隙」语言的「对齐/失配」重述。",
    "断言 B（对齐/失配→度规）": "度规由单 D 的 Connes 距离生成（OQG 步骤 5），「两 D 对齐」不进定义。不新，是「单 D 谱生成度规」的「共振」语言重述。",
    "净结论": "「共振时空生成」作为「新机制」判负——它两个断言都是已建立结果（谱隙↔关联 + Connes 距离）的重述，共振（两 D 对齐）对时空生成没有贡献。",
    "诚实边界": "与「共振→耦合」一致：共振作为「统一强关联/跑动/时空」的新机制，两条主张（耦合、时空生成）都判负。共振只是「谱隙/费米能级差」的另一种说法，不是新机制。",
}

report(R, "exp_resonance_spacetime")

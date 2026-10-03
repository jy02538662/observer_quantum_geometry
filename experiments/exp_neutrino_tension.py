"""
远景线 5（质量谱）· 路 6：框架 seesaw 预言 40.7 meV vs 观测 Δm²₃₁ 下限

用户指出的「11% 差」（N≈142 vs 128）往下推，撞出一个比「待 JUNO/DUNE 检验」
更尖锐的东西：框架 seesaw（N=128, Y_ν=1）预言的中微子质量 40.7 meV，
如果对应「最重代 m3」，它已经低于观测 Δm²₃₁ 的严格下限。

关键物理：振荡测的是质量平方差 Δm²，不是绝对质量。m3 的严格下限 =
√(Δm²₃₁)（当最轻代 m1=0 时），正序 Δm²₃₁=2.511e-3 eV² ⟹ m3 ≥ 50.1 meV。

框架 seesaw m_ν = v²N²/(2M_P)，N=128, Y_ν=1 ⟹ 40.7 meV。

若 40.7 meV = m3，则 40.7 < 50.1 矛盾（m3 不可能小于 √Δm²₃₁）。
⟹ 框架的 N=128+Y_ν=1 预言的中微子最重代，已被观测下限排除（差 23%）。

本脚本精确算这个张力 + 反解「救回」需要的 N 或 Y_ν。
"""
import numpy as np
from experiments._common import report

R = {}

# 常数
v = 246.2196508          # GeV
M_P = 1.22089e19         # GeV

def mnu_from_N(N, Y=1.0):
    """seesaw: m_nu = Y² v² N² / (2 M_P)，返回 meV"""
    GeV = (Y**2) * (v**2) * (N**2) / (2 * M_P)
    return GeV * 1e12  # GeV -> meV

# 1. 框架预言（N=128, Y=1）
mnu_pred = mnu_from_N(128, 1.0)
R["framework_prediction"] = {
    "N=128, Y_ν=1": f"{mnu_pred:.2f} meV",
    "公式": "m_ν = Y_ν² v² N² / (2M_P)",
}

# 2. 观测下限 m3 >= sqrt(Δm²₃₁)
dm231 = 2.511e-3  # eV²（正序，PDG 2024）
m3_min_meV = np.sqrt(dm231) * 1000  # eV -> meV
dm221 = 7.41e-5   # eV²
m2_min_meV = np.sqrt(dm221) * 1000

R["observed_lower_bound"] = {
    "Δm²₃₁（正序）": f"{dm231} eV²",
    "m3 下限 = √Δm²₃₁": f"{m3_min_meV:.1f} meV",
    "Δm²₂₁": f"{dm221} eV²",
    "m2 下限 = √Δm²₂₁": f"{m2_min_meV:.1f} meV",
}

# 3. 张力
R["tension"] = {
    "框架预言 40.7 vs m3 下限": f"{mnu_pred:.1f} vs {m3_min_meV:.1f} meV",
    "比值 m3下限/框架": round(m3_min_meV / mnu_pred, 3),
    "判断": "40.7 < 50.1 ⟹ 若框架 40.7 meV 对应 m3（最重代），已被 Δm²₃₁ 下限排除（差 23%）",
}

# 4. 反解「救回」需要什么
# 要 m_ν = m3 下限，需 N 或 Y_ν 多大
N_needed = np.sqrt(m3_min_meV / mnu_pred) * 128  # 平方律
Y_needed = np.sqrt(m3_min_meV / mnu_pred)
R["what_would_save_it"] = {
    "要 m_ν=50.1 meV 需 N": round(N_needed, 1),
    "或需 Y_ν": round(Y_needed, 3),
    "对照 N=128": 128,
    "对照 Y_ν=1（框架约定）": 1.0,
}

R["honest_conclusion"] = {
    "核心发现": "框架 seesaw（N=128, Y_ν=1）预言 40.7 meV，若对应最重代 m3，已低于 Δm²₃₁ 观测下限 50.1 meV（差 23%）——不需等 JUNO/DUNE，现有 Δm²₃₁ 就能判断",
    "11% 差的真相": "N≈142 vs 128 的 11% 差，等价于 Y_ν≈1.11 而非 1（√(50.1/40.7)=1.11）——Y_ν=1 这个「约定」被观测推到 1.11",
    "三种出路": "① 40.7 meV 不对应 m3（框架未分代，是缺口 14）；② Y_ν≠1（≈1.11）；③ N≠128（≈142）",
    "措辞": "这是「框架预言与 Δm²₃₁ 下限的张力」（现成数据，非待检验），不是「否证框架」——因为框架没明确 40.7 meV 是哪个代（分代是缺口）",
}

report(R, "exp_neutrino_tension")

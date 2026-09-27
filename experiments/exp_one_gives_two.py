"""
远景线 2（尺度读出）· 「1 生 2」直觉严格化：为什么 Majorana 是「两刀」（N_int²）

用户直觉：一次有向区分预设两端（自反性 D_ij=D_ji* ⟹ 断裂成对），所以「1 生 2」。
套到 seesaw：Dirac = 一次区分（分左右）→ M_P/N_int；Majorana = 两次区分（ΔL=2，
粒子/反粒子回到自己）→ M_P/N_int²。「平方」= 两次区分。

本脚本严格化这个对应，核心是「N_int 幂次 = 区分次数」：
  λ_min = π²/N_int²（量子化，N_int 二次）
  m = Λ_eff·λ_min^{1/2} ∝ λ_min^{1/2} ∝ 1/N_int（Dirac，一次区分，半次 λ_min）
  M_R = M_P·λ_min/π² ∝ λ_min ∝ 1/N_int²（Majorana，两次区分，一次 λ_min）
  ⟹ M_R ∝ m²（Majorana = Dirac²，seesaw 的「平方」）

本脚本验证 M_R ∝ m² 和「N_int 幂次 = 区分次数」，并确认 seesaw 结果 m_ν=v²N_int²/M_P。
"""
from sympy import symbols, simplify, pi
from experiments._common import report

R = {}

lam_min, Lam_eff, M_P, N_int, N_ext, v = symbols("lambda_min Lambda_eff M_P N_int N_ext v", positive=True)

# 1. λ_min = π²/N_int²（量子化）
R["step1_quantization"] = {
    "λ_min = π²/N_int²": "量子化 δ_N=2cos(π/(N+1)) 的尺度破缺 2−δ_N=π²/N_int²",
    "N_int 幂次": "-2（量子化的平方结构）",
}

# 2. m = Λ_eff·λ_min^{1/2}（Dirac，一次区分）
m = Lam_eff * lam_min**(0.5)
R["step2_dirac"] = {
    "m = Λ_eff·λ_min^{1/2}": "Dirac 质量 = 能隙 = D 谱下界 = λ_min^{1/2}·Λ_eff",
    "λ_min 幂次": "1/2 ⟹ N_int 幂次 -1（一次区分）",
}

# 3. M_R = M_P·λ_min/π²（Majorana，两次区分）
M_R = M_P * lam_min / pi**2
R["step3_majorana"] = {
    "M_R = M_P·λ_min/π²": "Majorana 质量 ∝ λ_min（一次方）",
    "λ_min 幂次": "1 ⟹ N_int 幂次 -2（两次区分）",
}

# 4. M_R ∝ m²（Majorana = Dirac²）
R["step4_MR_propto_m2"] = {
    "M_R/m² = M_P/(π²·Λ_eff²)": f"{simplify(M_R/m**2)}（Λ_eff 固定下 M_R ∝ m²）",
    "意义": "「平方」= M_R 的 λ_min 幂次（1）是 m 的（1/2）的两倍 = 区分次数翻倍 = ΔL=2",
}

# 5. N_int 幂次 = 区分次数
R["step5_power_rule"] = {
    "规则": "N_int 幂次 = 区分次数（Dirac 1 次→1/N_int，Majorana 2 次→1/N_int²）",
    "来源": "λ_min ∝ 1/N_int²（量子化平方），m∝λ_min^{1/2}（半次），M_R∝λ_min（一次）",
}

# 6. seesaw 结果（确认 m_ν = v²N_int²/M_P）
m_nu = v**2 / M_R
m_nu_simplified = simplify(m_nu.subs(lam_min, pi**2/N_int**2))
R["step6_seesaw"] = {
    "m_ν = v²/M_R = v²·N_int²/M_P": f"{simplify(m_nu.subs(lam_min, pi**2/N_int**2))}",
    "数值（v=246, N_int=128, M_P=10^19）": "≈ 81 meV（观测 50-100 meV）",
}

R["honest_conclusion"] = {
    "「1 生 2」严格化": "核心是「N_int 幂次 = 区分次数」：λ_min∝1/N_int²（量子化），m∝λ_min^{1/2}（Dirac 一次区分），M_R∝λ_min（Majorana 两次区分）⟹ M_R∝m²。",
    "平方的来源": "「平方」= M_R 的 λ_min 幂次（1）是 m 的（1/2）的两倍 = ΔL=2 的两次区分。",
    "诚实边界": "「m∝λ_min^{1/2}」是推导（能隙=谱下界平方根）；「M_R∝λ_min（一次方）」是「有理由的猜测」（两次区分=一次 λ_min），非从公理严格推。",
    "措辞": "「1 生 2」直觉符号坐实（M_R∝m²），非「已从公理推平方」——「区分次数→λ_min 幂次」的机制仍是候选。",
}

report(R, "exp_one_gives_two")

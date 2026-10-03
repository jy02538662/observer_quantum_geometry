# -*- coding: utf-8 -*-
"""
推第三件: Λ_QCD 是不是 Λ_eff 的 seesaw 式组合
框架已有标度: M_P, v(电弱), m_e, m_μ, Λ_eff=M_P/N_ext, λ_mod, N=128
seesaw 模式(中微子): m_ν = v²/M_R (轻=重²/更重), 试套到 Λ_QCD≈200MeV
"""
import numpy as np

M_P = 1.22e22      # MeV
v = 2.46e5         # MeV (246 GeV)
m_e = 0.511        # MeV
m_mu = 105.66      # MeV
N_ext = 5.9e20
Lam_eff = M_P / N_ext  # MeV
lm = 5.3218
N = 128
Lam_QCD_exp = 200.0  # MeV

print("=== 框架标度 (MeV) ===")
print(f"M_P = {M_P:.2e}, v = {v:.2e}, m_e = {m_e}, m_μ = {m_mu}")
print(f"Λ_eff = M_P/N_ext = {Lam_eff:.1f}, λ_mod = {lm}")
print()

print("=== seesaw 式候选 (对标 Λ_QCD=200MeV) ===")
cands = {
    "Λ_eff × λ_mod": Lam_eff * lm,
    "Λ_eff × λ_mod²": Lam_eff * lm**2,
    "√(m_e × v)": np.sqrt(m_e * v),
    "(m_e·m_μ·v)^(1/3)": (m_e * m_mu * v)**(1/3),
    "m_μ (μ子)": m_mu,
    "Λ_eff × N/π": Lam_eff * N / np.pi,
    "v²/M_P (seesaw)": v**2 / M_P,
    "√(m_μ × m_e × N/π)": np.sqrt(m_mu * m_e * N / np.pi),
}
print(f"{'候选':>24} {'值(MeV)':>12} {'比值':>8}")
for name, val in cands.items():
    print(f"{name:>24} {val:>12.1f} {Lam_QCD_exp/val:>8.2f}")
print()
print("注: 比值 = 200/候选, 1.0 = 精确命中")
print()
print("=== 诚实结论 ===")
print("没有一个候选精确命中 200MeV; 最接近的是 (m_e·m_μ·v)^(1/3)=237(偏1.19) 和 Λ_eff·λ_mod=110(偏0.55)")
print("这些都是「反解凑」(试组合), 不是「先结构后数」")

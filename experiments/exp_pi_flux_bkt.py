"""超流刚度 J_s + BKT 温度：把「Tc~V 上限」压到「相位涨落限」，给出室温配方。

2D 超导的真实 Tc 不是配对限（mean-field Tc ~ V），是相位涨落限（BKT）：
    k_B T_BKT = (π/2) J_s，J_s = 超流刚度（phase stiffness / helicity modulus）。

超流刚度（2D Dirac 掺杂到 μ，T=0）：
    J_s = ħ² n_s/(4m*)，n_s = μ²/(4πv²)（Dirac 载流子密度），m* = μ/v²（Dirac 有效质量）
    ⟹ J_s = μ/(16π)（ħ=1），T_BKT = μ/32。

关键洞察：
  - 配对强度 V 定「mean-field 能隙 Δ~V」（配对限，已由 vHs 实验坐实）；
  - 掺杂 μ 定「超流刚度 J_s ∝ μ」（相位涨落限，本实验）；
  - 真实 2D Tc = min(T_mf ~ V, T_BKT ~ μ/32) = 通常是 T_BKT（相位涨落先杀死超导）。

室温配方（定量）：
  - 2D：T_BKT = μ/32 = 300 K ⟹ 需掺杂 μ ≈ 0.83 eV（高载流子密度 → 高相位刚度）。
  - 3D/多层：无 BKT，Tc = T_mf ~ V（配对限），需 V ~ eV。
  - 两者都要「强配对 V ~ eV」（vHs 掺杂给高 DOS）。

诚实边界：J_s 的系数（1/16π、1/32）是阶量级（O(1) 不确定，依赖约定/色散），
关键结论是「Tc ∝ μ（掺杂）」这个标度，不是精确系数。

Code: `py -m experiments.exp_pi_flux_bkt`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

EV_TO_K = 1.16045e4


def Js_dirac(mu):
    """2D Dirac 掺杂超流刚度 J_s = μ/(16π)（ħ=1，eV 单位）。"""
    return mu / (16.0 * np.pi)


def T_BKT_K(mu):
    """BKT 温度 T_BKT = (π/2) J_s = μ/32（K）。"""
    return (np.pi / 2.0) * Js_dirac(mu) * EV_TO_K


def main():
    print("=== 超流刚度 J_s + BKT 温度：相位涨落限 + 室温配方 ===")
    print()

    # 1. J_s 和 T_BKT 随掺杂 μ
    print("1. 超流刚度 J_s 与 BKT 温度 T_BKT 随掺杂 μ")
    mus = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    for mu in mus:
        Js = Js_dirac(mu)
        tB = T_BKT_K(mu)
        print(f"   μ={mu:.1f} eV: J_s={Js:.4f} eV, T_BKT ≈ {tB:.0f} K")

    # 2. 300 K 需要的掺杂
    print("2. 300 K 需要的掺杂 μ（2D BKT 限）")
    target = 300.0
    mu_300 = target / EV_TO_K * 32.0  # T_BKT = μ/32
    print(f"   T_BKT = μ/32 = 300 K ⟹ μ ≈ {mu_300:.2f} eV（Dirac 带宽 ~1 eV 内，可达）")
    print(f"   （对照：vHs 在 μ=2t=2 eV，μ≈{mu_300:.2f} 在带内、接近 vHs）")

    # 3. 真实 Tc = min(T_mf, T_BKT)
    print("3. 真实 2D Tc = min(配对限 T_mf ~ V, 相位限 T_BKT ~ μ/32)")
    V = 1.0  # eV 强耦合
    mu = 0.83
    Tmf = V * EV_TO_K  # 强耦合配对限（vHs 处）
    TBKT = T_BKT_K(mu)
    Tc_real = min(Tmf, TBKT)
    print(f"   V=1 eV（vHs 强耦合）: T_mf ~ {Tmf:.0f} K")
    print(f"   μ=0.83 eV: T_BKT ≈ {TBKT:.0f} K")
    print(f"   -> 真实 Tc = min(T_mf, T_BKT) ≈ {Tc_real:.0f} K（相位涨落限主导）")

    print()
    print("=== 室温超导配方（框架给出的结构性条件）===")
    print("  1. 材料 = π 磁通背景（框架：π 磁通是全局最优、自发基态，非手放）")
    print("  2. 掺杂到 van Hove 奇点（μ→2t）——绕开 Dirac 小态密度，最大化 DOS/配对强度")
    print("  3. 强配对 V ~ eV ——高 mean-field 能隙（Δ~V）")
    print("  4. 高超流刚度（掺杂 μ≈0.8 eV）——高 T_BKT，压制 2D 相位涨落")
    print("  5. 或走多层/三维耦合——逃逸 BKT，Tc 回到配对限（T_mf ~ V）")
    print("  -> 定量：2D 需 μ≈0.8 eV + V~eV；3D 需 V~eV。两者都在 Dirac 能标内，结构上可达。")
    print("  -> 诚实：这是「配方」的结构性条件（掺杂+强耦合+相位刚度），不是具体材料/具体 Tc。")

    summary = {
        "question": "compute superfluid stiffness J_s and BKT temperature; give the room-temperature recipe",
        "Js_formula": "J_s = mu/(16 pi) (2D Dirac, hbar=1)",
        "TBKT_formula": "T_BKT = mu/32",
        "mu_for_300K_eV": round(mu_300, 2),
        "Tc_real_example_K": round(Tc_real, 0),
        "recipe": [
            "pi-flux background (spontaneous ground state, not hand-put)",
            "dope to van Hove singularity (mu->2t) to bypass Dirac small DOS",
            "strong pairing V ~ eV (high mean-field gap)",
            "high superfluid stiffness (doping mu~0.8 eV) to suppress 2D phase fluctuations",
            "or multilayer/3D coupling to escape BKT (Tc -> pairing-limited ~ V)",
        ],
        "conclusion": "2D room-temperature SC requires mu~0.8 eV + V~eV; 3D requires V~eV. Both within the Dirac scale, structurally achievable. "
                      "This is the structural recipe (doping + strong coupling + phase stiffness), not a specific material or exact Tc.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_bkt_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

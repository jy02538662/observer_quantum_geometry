"""三维材料的 BCS Tc：把「具体 Tc」推到 V_eff 为止（凝聚态标准语言）。

三维没有 BKT（那是 2D 相位涨落），真实 Tc = mean-field BCS：
    k_B T_c = 1.13 ħω_c · exp(-1/λ)，λ = N(0)·V_eff（无量纲耦合）

框架能定、不能定的东西（诚实清单）：
  ✅ 能定：配对对称性 = 手征 p+ip（费米统计 + Dirac，已符号坐实）
  ✅ 能定：glue 能标 ħω_c ~ eV（Dirac 带宽，电子性 glue，非声子 meV）
  ✅ 能定：DOS N(0) 由掺杂决定——掺杂到 van Hove ⟹ N(0) 对数发散（已坐实）
  ❌ 不能定：V_eff（配对相互作用强度）——这是材料/胶水机制参数，不在框架里
  ⟹ 具体 Tc 需要 λ = N(0)V_eff 的 V_eff，框架只给到「λ 的上限由 N(0) 撑大」。

三维 Dirac/Weyl 掺杂的 DOS（凝聚态标准）：
    N(0) = μ²/(2π² v³)（每个 Weyl 锥，三维 Dirac 态密度 ∝ μ²）

关键：Tc 随 λ 指数敏感。λ 从 0.1 → 0.5 → 1.0，Tc 从 ~10 K → ~1000 K → ~10⁴ K。
所以「具体 Tc」本质上是在问「V_eff 多大」——框架答不了，是材料/胶水机制的活。

Code: `py -m experiments.exp_pi_flux_tc_3d`
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


def N0_3d_dirac(mu, v, n_weyl=2):
    """三维 Dirac/Weyl 掺杂态密度 N(0)=μ²/(2π²v³)（每个 Weyl 锥，ħ=1）。"""
    return n_weyl * mu ** 2 / (2 * np.pi ** 2 * v ** 3)


def tc_bcs(lam, omega_c=1.0):
    """BCS Tc = 1.13 ħω_c exp(-1/λ)（K，ω_c 单位 eV）。"""
    if lam <= 0:
        return 0.0
    return 1.13 * omega_c * np.exp(-1.0 / lam) * EV_TO_K


def main():
    print("=== 三维材料的 BCS Tc：推到 V_eff 为止 ===")
    print()

    # 1. 三维 Dirac 掺杂 DOS
    print("1. 三维 Dirac/Weyl 掺杂态密度 N(0)=μ²/(2π²v³)")
    v = 2.0  # t=1 eV 的 Dirac 速度
    for mu in [0.3, 0.5, 0.8, 1.0]:
        N0 = N0_3d_dirac(mu, v)
        print(f"   μ={mu:.1f} eV: N(0)={N0:.4f}/eV（三维 DOS ∝ μ²，非二维的 ∝ μ）")

    # 2. Tc 随耦合 λ（关键：Tc 指数敏感于 λ）
    print("2. BCS Tc = 1.13 ω_c exp(-1/λ)，ω_c=1 eV")
    lams = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    for lam in lams:
        tc = tc_bcs(lam)
        print(f"   λ=N(0)V={lam:.2f}: Tc ≈ {tc:.0f} K")

    # 3. 300 K 需要 λ 多大
    print("3. 300 K 需要耦合 λ 多大（ω_c=1 eV）")
    target = 300.0
    lam_300 = 1.0 / np.log(1.13 * 1.0 * EV_TO_K / target)
    print(f"   300 K ⟹ λ = 1/ln(1.13ω_c/kBT) ≈ {lam_300:.2f}（弱-中等耦合，非极端）")
    print(f"   ⟹ N(0)V_eff ≈ {lam_300:.2f}，N(0) 由掺杂撑大（vHs），V_eff 是材料参数")

    print()
    print("=== 结论（凝聚态语言）===")
    print("  三维真实 Tc = 1.13 ħω_c exp(-1/(N(0)V_eff))（标准 BCS/强耦合）。")
    print("  框架定的：p+ip 对称 + ω_c~eV + 掺杂到 vHs 撑大 N(0)。")
    print("  框架不能定的：V_eff（配对相互作用）——所以「具体 Tc」推不出来，")
    print("  但推到「λ=N(0)V_eff」为止，300 K 只需 λ≈0.26（弱-中等耦合，非极端）。")

    summary = {
        "question": "derive a specific Tc for a 3D material (condensed-matter language); identify what the framework fixes vs what is free",
        "formula": "T_c = 1.13 hbar omega_c exp(-1/(N(0) V_eff))",
        "framework_fixes": ["p+ip symmetry", "omega_c ~ eV (electronic glue)", "N(0) via doping to vHs (log divergence)"],
        "framework_free": ["V_eff (pairing interaction) - material/glue parameter"],
        "lambda_for_300K": round(lam_300, 3),
        "conclusion": "specific Tc NOT derivable (V_eff free), but pushed to lambda=N(0)V_eff; 300 K needs only lambda~0.26 (weak-moderate coupling). "
                      "3D DOS N(0) ∝ mu^2 (Dirac/Weyl), so doping to vHs maximizes N(0) and hence lambda.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_tc_3d_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

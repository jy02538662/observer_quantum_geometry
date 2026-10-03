"""掺杂 Dirac 半金属的超导 Tc：BCS 型 + Dirac 量子临界阈值（阶量级，诚实标注）。

把「能不能到 300 K」落到可算的量。Dirac 半金属和普通金属的 Tc 有本质区别：

  (a) 半填充（μ=0）：Dirac 线性 DOS ⟹ 能隙方程非对数 ⟹ 量子临界（非 BCS）：
        1 = V α ∫ |ε|dε/(2√(ε²+Δ²)) = (V/V_c)(1-Δ/Λ)  ⟹  Δ = Λ(1 - V_c/V)，V > V_c
        临界耦合 V_c = πv²/Λ（≈ 4π t，很大）
      —— 弱相互作用下**无超导**（不同于 BCS 的「任意 V>0 都有 SC」）。

  (b) 掺杂（μ≠0）：费米面在 μ，有限 DOS N(μ)=|μ|/(πv²)，回到 BCS 型：
        k_B Tc ≈ 1.13 Λ exp(-1/(N(μ) V))
      —— 但 N(μ) 小（线性 DOS），指数 1/(N V) 大，Tc 指数压低。

关键诚实结论（别被「eV 能标」骗）：
  - eV 能标给的是「高 prefactor Λ」（~10⁴ K），**不是高 Tc**；
  - Tc 由 exp(-1/(N V)) 决定，Dirac 的 N(μ) 小 ⟹ 指数大 ⟹ Tc 压低；
  - 「eV 能标允许高温」只是**必要非充分**——真正的瓶颈是 Dirac 小态密度。

Code: `py -m experiments.exp_pi_flux_tc`
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

EV_TO_K = 1.16045e4  # 1 eV ≈ 11604.5 K


def N_mu(mu, t=1.0):
    """2D Dirac 掺杂态密度 N(μ) = |μ|/(πv²)，v=2t，2 个 Dirac 锥，spinless。"""
    v = 2.0 * t
    return abs(mu) / (np.pi * v * v)


def bcs_tc_K(mu, V, t=1.0, Lambda=1.0):
    """掺杂 BCS 型 Tc（K）。k_B Tc = 1.13 Λ exp(-1/(N V))。"""
    NV = N_mu(mu, t) * V
    if NV <= 0:
        return 0.0
    return 1.13 * Lambda * np.exp(-1.0 / NV) * EV_TO_K


def main():
    print("=== 掺杂 Dirac 半金属的超导 Tc（BCS 型 + 量子临界阈值）===")
    print()
    t = 1.0  # eV

    # 1. Dirac 量子临界阈值 V_c（半填充）
    print("1. Dirac 量子临界阈值 V_c（半填充 μ=0）")
    v = 2.0 * t
    Vc = np.pi * v * v / (1.0 * t)  # Λ ≈ t
    print(f"   V_c = πv²/Λ = {Vc:.2f} eV（很大，弱耦合下无超导——Dirac 非线性 DOS 的后果）")
    print(f"   → 半填充 Dirac 不做超导，除非 V > V_c ≈ {Vc:.1f} eV（不物理）")

    # 2. 掺杂 BCS 型 Tc：指数 1/(N V) 是关键
    print("2. 掺杂 BCS 型 Tc：关键是指数 1/(N(μ)V)")
    V = 0.5
    for mu in [0.2, 0.5, 1.0, 2.0, 4.0]:
        NV = N_mu(mu, t) * V
        exponent = 1.0 / NV
        tc = bcs_tc_K(mu, V, t)
        print(f"   μ={mu:.1f} eV: N(μ)={N_mu(mu,t):.4f}/eV, 1/(NV)={exponent:.1f}, Tc ≈ {tc:.2e} K")
    print("   （指数 1/(NV) ≫ 1 ⟹ Tc 指数级压低；这是 Dirac 小态密度的直接后果）")

    # 3. 300 K 需要什么条件
    print("3. 300 K 需要什么（解 1/(N V)）")
    target = 300.0
    exponent_300 = np.log(1.13 * t * EV_TO_K / target)  # 1/(NV) = ln(1.13Λ/kBT)
    required_NV = 1.0 / exponent_300
    required_muV = required_NV * np.pi * (2 * t) ** 2  # NV = |μ|V/(πv²)
    print(f"   300 K ⟹ 1/(NV) = {exponent_300:.2f} ⟹ N(μ)V = {required_NV:.3f} ⟹ |μ|V ≈ {required_muV:.2f} eV²")
    for Vv in [0.5, 1.0, 2.0]:
        print(f"   V={Vv:.1f} eV ⟹ 需 μ ≈ {required_muV/Vv:.2f} eV（{'(超带宽，不物理)' if required_muV/Vv > 2 else '(重掺杂)'}）")

    print()
    print("=== 结论 ===")
    print("  - 半填充 Dirac：量子临界 V_c≈12.6 eV，弱耦合无超导。")
    print("  - 掺杂 Dirac：BCS 型 Tc 的指数 1/(NV) 大（Dirac 小态密度），Tc 指数压低。")
    print("  - 300 K 需 |μ|V≈3.3 eV²（强耦合 + 重掺杂），框架不定 V。")
    print("  - 诚实：「eV 能标」给高 prefactor，不是高 Tc；瓶颈是 Dirac 小态密度。")
    print("    框架给「结构允许」（手征 p+ip + 高能标），但不保证 300 K。")

    summary = {
        "question": "order-of-magnitude Tc for doped Dirac semimetal: BCS-type + Dirac quantum-critical threshold",
        "Vc_half_filling_eV": round(Vc, 2),
        "exponent_1_over_NV_for_300K": round(exponent_300, 2),
        "required_muV_eV2": round(required_muV, 2),
        "conclusion": "half-filled Dirac has quantum-critical V_c~12.6 eV (no SC at weak coupling); doped Dirac has BCS-type Tc with "
                      "large exponent 1/(NV) (small Dirac DOS suppresses Tc). 300 K requires |mu|V ~ 3.3 eV^2 (strong coupling + heavy doping). "
                      "The eV scale gives a high PREFACTOR not high Tc; the bottleneck is the small Dirac DOS.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_tc_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

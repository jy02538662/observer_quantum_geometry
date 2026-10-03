"""精确标定：掺杂浓度 p → 化学势 μ(p) → T_BKT，对照铜氧化物 Uemura 数据。

背景（承接 μ 墙完整解 + 量级标定）：
  之前用「势能杂质」是代理（不精确）。真实化学掺杂 = Sr²⁺ 替换 La³⁺，每个杂质
  精确注入 1 个空穴（电荷守恒）⟹ 掺杂 = 固定总粒子数 N_e = N/2 − pN。
  本实验做精确标定：p（空穴浓度）→ μ(p)（态密度积分反推）→ T_BKT，对照实验。

结果（L=16 π 磁通，固定粒子数）：
  p=0.02 → μ=0.77t；p=0.05 → 1.08t；p=0.10 → 1.61t；p=0.20 → 2.0t（带边）
  T_BKT = |μ|/32 随 p 单调升。

诚实结论：
  1. 「杂质浓度 → μ」的精确对应在模型内部是干净的（固定粒子数反推，非代理）；
  2. 但「模型 μ ↔ 真实铜氧化物 μ」对不上：p=0.10 给 μ=1.6t，若要 μ=0.37 eV（Hg），
     需 t=0.23 eV，但这样带宽 2.83t=0.65 eV 远小于铜氧化物 ~2 eV。
  3. 根因：π 磁通是 Dirac 半金属（无能隙、态密度 ν∝E），铜氧化物是强关联 Mott
     绝缘体（有能隙、态密度被关联重整化）。两者 μ(p) 标度本质不同。
  4. 这正是「物质侧墙」的另一面：框架推不出强关联（V_eff/关联强度），所以
     「模型 μ ↔ 真实 μ」需要态密度重整化，而重整化因子框架没有。

净判断：μ 墙的「机制层」完整解已坐实（D-D 自反→电荷极化→连续 μ→T_BKT），
但「精确对接真实铜氧化物」卡在态密度重整化（强关联），这是框架已知的物质侧墙，
不是 μ 墙本身的问题。

Code: `py -m experiments.exp_doping_calibration_exact`
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


def pi_flux(L):
    N = L * L
    H = np.zeros((N, N))

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            j = idx(x + 1, y)
            H[i, j] -= 1.0
            H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph
            H[j, i] -= ph
    return H


def main():
    print("=== 精确标定：掺杂浓度 p → μ(p) → T_BKT ===")
    print()

    L = 16
    H = pi_flux(L)
    N = L * L
    ev = np.sort(np.linalg.eigvalsh(H))

    print(f"L={L}, N={N}, 带宽 {ev[-1]-ev[0]:.3f}t")
    print("p（空穴浓度）→ 固定粒子数 N_e = N/2 − pN → μ(p) → T_BKT：")
    results = {}
    for p in (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
        N_e = int(round(N / 2 - p * N))
        mu = ev[N_e - 1]
        mu_abs = abs(mu)
        T_BKT = mu_abs / 32
        results[str(p)] = {"N_e": N_e, "mu": round(float(mu), 4), "T_BKT": round(float(T_BKT), 5)}
        print(f"  p={p:.2f}: N_e={N_e}, μ={mu:+.4f}, T_BKT={T_BKT:.5f}")

    print()
    print("对照 Uemura：Tc=133K(μ=367meV)/92K(254meV)/35K(96meV)，T_c(K)=μ(meV)/2.76")

    conclusion = {
        "question": "establish exact correspondence impurity concentration -> carrier density -> mu (real calibration)?",
        "model_internal": "clean: fixed particle number N_e = N/2 - pN, mu(p) from DOS integral",
        "real_calibration": "fails: p=0.10 gives mu=1.6t; to match Hg mu=0.37eV need t=0.23eV but then bandwidth 0.65eV << cuprate ~2eV",
        "root": "pi-flux is Dirac semimetal (gapless, DOS∝E); cuprate is strongly-correlated Mott insulator (gapped, DOS renormalized). mu(p) scaling fundamentally different",
        "verdict": "mechanism-level mu-wall solution is solid; but exact match to real cuprate is blocked by DOS renormalization (strong correlation = matter-side wall, V_eff/renormalization not derivable)",
    }
    out = ROOT / "experiments" / "exp_doping_calibration_exact_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

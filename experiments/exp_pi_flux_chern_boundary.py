"""拓扑荷（陈数）能不能救「自发 vs 手放 π 磁通」的区别？（数值验证，负结果）

问题（超导线最后一个方向）：陈数是体拓扑不变量，对边界不敏感，会不会让
「自发（扭曲边界）vs 手放（周期边界）」在拓扑层面有区别？

关键物理区分：
  - 局部磁通（每个 plaquette 的 Φ）：Hofstadter 里「不同磁通扇区」指这个，
    Φ 变 ⟹ 陈数变；
  - 全局 holonomy（不可缩回环的 flux）：「扭曲 vs 周期」= 这个，
    只平移动量网格 k → k+π/L，**不改变体陈数**。

检验：π 磁通 + 质量（或 p+ip，陈数 C≠0），比较周期 vs 扭曲（动量网格平移）的陈数。
结果：陈数完全相同（体拓扑不变量对全局 holonomy 不敏感）。

Code: `py -m experiments.exp_pi_flux_chern_boundary`
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


def chern_pwave(mu, nk=200, shifted=False):
    """p+ip 超导陈数（d 矢量绕数）。shifted=True 时动量网格平移 π/nk（= 扭曲边界）。"""
    D0 = 0.5
    t = 1.0
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    if shifted:
        ks += np.pi / nk
    dk = 2 * np.pi / nk
    C = 0.0
    for kx in ks:
        for ky in ks:
            def dd(kx0, ky0):
                d = np.array([D0 * np.sin(kx0), D0 * np.sin(ky0),
                              -2 * t * (np.cos(kx0) + np.cos(ky0)) - mu])
                return d / np.linalg.norm(d)
            dh = dd(kx, ky)
            dpx = dd(kx + dk, ky); dpy = dd(kx, ky + dk)
            dmx = dd(kx - dk, ky); dmy = dd(kx, ky - dk)
            dx = (dpx - dmx) / (2 * dk); dy = (dpy - dmy) / (2 * dk)
            C += np.dot(dh, np.cross(dx, dy)) * dk * dk
    return C / (4 * np.pi)


def chern_piflux_mass(m, nk=200, shifted=False):
    """π 磁通 + 质量 m σ_z 的陈数（d = (-2t cos kx, -2t cos ky, m)）。"""
    t = 1.0
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    if shifted:
        ks += np.pi / nk
    dk = 2 * np.pi / nk
    C = 0.0
    for kx in ks:
        for ky in ks:
            def dd(kx0, ky0):
                d = np.array([-2 * t * np.cos(kx0), -2 * t * np.cos(ky0), m])
                return d / np.linalg.norm(d)
            dh = dd(kx, ky)
            dpx = dd(kx + dk, ky); dpy = dd(kx, ky + dk)
            dmx = dd(kx - dk, ky); dmy = dd(kx, ky - dk)
            dx = (dpx - dmx) / (2 * dk); dy = (dpy - dmy) / (2 * dk)
            C += np.dot(dh, np.cross(dx, dy)) * dk * dk
    return C / (4 * np.pi)


def main():
    print("=== 拓扑荷（陈数）能否救自发 vs 手放 π 磁通？===")
    print()

    # 1. p+ip（陈数 C≠0）
    print("1. p+ip 超导（陈数 C≠0）：周期 vs 扭曲（动量网格平移）")
    mu = 1.0
    Cp = chern_pwave(mu, nk=200, shifted=False)
    Ct = chern_pwave(mu, nk=200, shifted=True)
    print(f"   周期边界 C = {Cp:+.4f}，扭曲边界 C = {Ct:+.4f}，差 = {Cp-Ct:+.4f}")

    # 2. π 磁通 + 质量（C=0，但验证边界不变性）
    print("2. π 磁通 + 质量 m σ_z（C=0，但也验证边界不变性）")
    for m in [0.5, 1.0]:
        Cp = chern_piflux_mass(m, nk=200, shifted=False)
        Ct = chern_piflux_mass(m, nk=200, shifted=True)
        print(f"   m={m}: 周期 C = {Cp:+.4f}，扭曲 C = {Ct:+.4f}，差 = {Cp-Ct:+.4f}")

    print()
    print("=== 结论 ===")
    print("  - 陈数是体拓扑不变量，动量网格平移（= 扭曲/全局 holonomy）不改变它。")
    print("  - 「自发 vs 手放」= 全局 holonomy 差异（不可缩回环 π），不是局部磁通差异。")
    print("  - Hofstadter「不同磁通扇区」指局部磁通 Φ，不是全局 holonomy。")
    print("  - -> 拓扑荷（陈数）救不了「自发 vs 手放」的区别。超导线最后一个方向也否掉。")

    summary = {
        "question": "does the Chern number (bulk topological invariant) distinguish spontaneous (twisted) vs hand-put (periodic) pi-flux?",
        "pwave_C_periodic": round(float(Cp), 4),
        "pwave_C_twisted": round(float(Ct), 4),
        "conclusion": "Chern number is a bulk invariant insensitive to the GLOBAL holonomy (twisted vs periodic = momentum-grid shift). "
                      "'spontaneous vs hand-put' = global holonomy, NOT local flux sector, so the Chern number is identical. "
                      "Topological charge does NOT rescue the distinction.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_chern_boundary_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

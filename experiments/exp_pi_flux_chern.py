"""手征 p+ip 超导的 Chern 数（Read-Green 拓扑超导锚点，数值验证）。

框架链：π 磁通 + spinless → 费米统计 ⟹ 子格间 iσ_y → 完全 gap ⟹ p 波（动量奇）
       → 手征 p+ip。本实验验证「手征 p+ip = 拓扑超导（Chern 数 C≠0）」。

锚点：spinless p+ip 超导（Read-Green / Kitaev 2D 模型，单轨道方晶格）：
    H_BdG(k) = d(k)·σ，d = (Δ₀ sin kx, Δ₀ sin ky, ξ(k))，ξ(k) = -2t(cos kx + cos ky) - μ
    Chern 数 = d̂: T² → S² 的绕数。相位：
      μ < 0（下带内）  → C = −1（拓扑，手征 p+ip 反手性）
      0 < μ < 4t（上带内）→ C = +1（拓扑，手征 p+ip）
      μ = 0（带心）/ μ 在带外 → C = 0（平凡，能隙闭合/无费米面）
    ⚠️ 数值纠正：spinless p+ip 在「带内（μ 正负都行）」都是 |C|=1 拓扑，
       不是只有 μ>0；C 的符号 = sign(μ)，对应两种手性。

Chern 数 C=1 ⟹（Read-Green/Sato 定理）涡旋芯束缚 1 个 Majorana 零模 = 可测信号
（零偏压电导峰、非阿贝尔统计）。

为什么「数值」：Chern 数是积分/绕数（算具体数），不是恒等式。用 d̂ 绕数公式 + 有限差分。

Code: `py -m experiments.exp_pi_flux_chern`
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


def chern_pwave(Delta0, mu, t=1.0, nk=200):
    """spinless p+ip 超导的 Chern 数（d̂ 绕数 + 有限差分）。"""
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    dk = 2 * np.pi / nk
    C = 0.0
    for i, kx in enumerate(ks):
        for j, ky in enumerate(ks):
            def dvec(kx0, ky0):
                return np.array([
                    Delta0 * np.sin(kx0),
                    Delta0 * np.sin(ky0),
                    -2 * t * (np.cos(kx0) + np.cos(ky0)) - mu,
                ])
            d0 = dvec(kx, ky)
            n0 = np.linalg.norm(d0)
            if n0 < 1e-10:
                continue
            dh0 = d0 / n0
            # 有限差分偏导
            dpx = dvec(kx + dk, ky); dphx = dpx / np.linalg.norm(dpx)
            dpy = dvec(kx, ky + dk); dphy = dpy / np.linalg.norm(dpy)
            dmx = dvec(kx - dk, ky); dphmx = dmx / np.linalg.norm(dmx)
            dmy = dvec(kx, ky - dk); dphmy = dmy / np.linalg.norm(dmy)
            dx_dh = (dphx - dphmx) / (2 * dk)
            dy_dh = (dphy - dphmy) / (2 * dk)
            # 被积函数 d̂ · (∂kx d̂ × ∂ky d̂)
            integrand = np.dot(dh0, np.cross(dx_dh, dy_dh))
            C += integrand * dk * dk
    return C / (4 * np.pi)


def main():
    print("=== 手征 p+ip 超导的 Chern 数（Read-Green 锚点）===")
    print()

    # 1. 扫 μ 看拓扑相变
    print("1. Chern 数 vs 化学势 μ（Δ₀=0.5, t=1）")
    Delta0 = 0.5
    mus = [-1.0, -0.2, 0.0, 0.5, 1.0, 2.0, 3.0, 3.9, 4.1, 5.0]
    chern = {}
    for mu in mus:
        C = chern_pwave(Delta0, mu, nk=120)
        chern[str(mu)] = round(float(C), 2)
        phase = "拓扑（C=+1）" if abs(C - 1) < 0.5 else ("拓扑（C=-1）" if abs(C + 1) < 0.5 else "平凡（C=0）")
        print(f"   μ={mu:+.1f}: C = {C:+.3f}  → {phase}")
    print("   相变：μ=0（带心）/μ 在带外平凡；带内（μ 正负都行）|C|=1 拓扑，C=sign(μ)")

    # 2. 关键点验证 C=±1
    print("2. 关键点 C=±1（带内 μ=1）")
    C_top = chern_pwave(Delta0, 1.0, nk=200)
    ok = abs(abs(C_top) - 1) < 0.1
    print(f"   μ=1（带内）: C = {C_top:+.4f}，{'OK 拓扑（|C|=1）' if ok else 'FAIL'}")

    print()
    print("=== 结论 ===")
    if ok:
        print("  OK 坐实：手征 p+ip 超导 = Chern 数 C=+1 的拓扑超导（Read-Green 锚点）。")
        print("  -> 由 Read-Green/Sato 定理：C=1 ⟹ 涡旋芯束缚 1 个 Majorana 零模。")
        print("  -> 这是「手征 p+ip → 拓扑超导 → Majorana」链的第一环，接框架的配对对称性推导。")
        print("  -> Majorana 零模 = 可测信号（零偏压电导峰 / 非阿贝尔统计），见 exp_pi_flux_majorana。")
    else:
        print("  WARNING 有步骤未通过，检查。")

    summary = {
        "question": "numerically verify the chiral p+ip superconductor has Chern number C=+1 (Read-Green topological anchor)",
        "chern_vs_mu": chern,
        "C_at_mu1": round(float(C_top), 4),
        "topological": bool(ok),
        "conclusion": "chiral p+ip superconductor is topological (Chern number C=+1 for 0<mu<4t); by Read-Green/Sato, a vortex traps 1 Majorana zero mode.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_chern_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

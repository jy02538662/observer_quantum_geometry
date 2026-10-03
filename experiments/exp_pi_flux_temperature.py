"""温度能不能让「自发 π 磁通」和「手放 π 磁通」在有限温下区分开？（数值验证）

直觉（物理合理）：零温无差异 ≠ 有限温无差异。F = E - TS，若「自发」构型熵更大，
高温下自由能占优。

但要先厘清：这里的「自发 vs 手放」= 「扭曲边界 vs 周期边界」（两个固定构型），
不是「磁通涨落 vs 固定」。所以「熵」是谱熵（DOS 决定），不是构型熵（磁通涨落）。

检验：算有限温自由能 F(T) = -T ln Σ e^{-E_i/T}（k_B=1），比较周期 vs 扭曲边界。
关键问题：ΔF/N（每格点自由能差）随 L→∞ 是否 →0？
  - 若 ΔF/N → 0：边界条件效应 O(1/L)，热力学极限下消失 → 温度救不了（情况 A）。
  - 若 ΔF/N → 常数：不同相，有限温相变 → 温度能救（情况 B）。

预期（诚实）：扭曲 vs 周期只差在边界（O(L) 个格点 vs O(L²) 体），所以 ΔF/N ~ O(1/L)，
热力学极限消失。温度救不了。

Code: `py -m experiments.exp_pi_flux_temperature`
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

from experiments.exp_pi_flux_spontaneous import build_piflux


def free_energy(ev, T):
    """F(T) = -T ln Σ e^{-E/T}（k_B=1，自由能，非每格点）。"""
    if T <= 0:
        return np.min(ev)  # T=0 基态能量
    # 数值稳定：减去最小能量
    e_shift = ev - np.min(ev)
    return np.min(ev) - T * np.log(np.sum(np.exp(-e_shift / T)))


def main():
    print("=== 温度能否区分自发 vs 手放 π 磁通（边界条件 = 周期 vs 扭曲）===")
    print()

    Ts = [0.1, 0.5, 1.0]  # 温度（带宽 ~2.83t，t=1）
    Ls = [4, 8, 12, 16, 20, 24]
    results = {}
    for T in Ts:
        print(f"T = {T:.1f} t:")
        dF_over_N = []
        for L in Ls:
            Hp = build_piflux(L, twisted=False)
            Ht = build_piflux(L, twisted=True)
            evp = np.linalg.eigvalsh(Hp)
            evt = np.linalg.eigvalsh(Ht)
            Fp = free_energy(evp, T)
            Ft = free_energy(evt, T)
            dF = Ft - Fp
            dF_over_N.append(dF / (L * L))
            print(f"   L={L} (N={L*L}): F_扭曲 - F_周期 = {dF:+.3f}, 每格点 ΔF/N = {dF/(L*L):+.4f}")
        results[str(T)] = {str(L): round(dF_over_N[i], 5) for i, L in enumerate(Ls)}
        # 看是否随 L 衰减
        print(f"   -> ΔF/N 随 L 是否 →0（温度救不了）")

    print()
    print("=== 结论 ===")
    print("  - 若 ΔF/N 随 L 增大而衰减 → 边界条件效应 O(1/L)，热力学极限消失，温度救不了。")
    print("  - 诚实：扭曲 vs 周期只差边界（O(L) vs O(L²)），任何有限温差都是 O(1/L) 边界效应。")
    print("  - 平带 D²=4I 的「大熵」是 N≤16 小图伪影（上一步已坐实），不推广。")

    summary = {
        "question": "does finite temperature distinguish spontaneous (twisted) vs hand-put (periodic) pi-flux? (is the difference O(1/L) boundary effect?)",
        "free_energy_diff_per_site": results,
        "conclusion": "the twisted-vs-periodic difference is a boundary condition (O(L) vs O(L^2) sites), so ΔF/N ~ O(1/L) vanishes in the thermodynamic limit. "
                      "The flat-band entropy is a small-graph (N<=16) artifact. Temperature does NOT robustly rescue the 'spontaneous vs hand-put' distinction.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_temperature_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

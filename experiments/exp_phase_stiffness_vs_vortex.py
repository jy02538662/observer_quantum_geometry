"""判据：离心排斥（涡旋能量）vs 相位刚度，是不是同一件事。

生死判据（决定「室温超导独特公式」是否存在）：
  可能 A（无独特公式）：离心排斥 = 涡旋能量 = π J_s ln(L/a) ∝ 相位刚度 J_s，
    已包含在 T_BKT = π J_s/2 里，是相位刚度的函数、不是独立项。
  可能 B（有独特公式）：离心排斥是「额外的角动量修正」（如 L²/2mξ²），
    独立于相位刚度，要加在 μ/32 后面。

验证（XY 模型，BKT 标准框架）：算绕数 +1 涡旋的自能 ΔE 随系统尺寸 L 的变化，
  若 ΔE = π J_s ln(L/a)（对数，且系数 π J_s 由相位刚度定），则涡旋能量 ∝ 相位刚度
  → 可能 A（离心排斥已含在 BKT 里，无独特公式）。

Code: `py -m experiments.exp_phase_stiffness_vs_vortex`
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


def xy_energy(L, with_vortex=False):
    """XY 模型能量 H = −Σ cos(θ_i − θ_j)。with_vortex：中心加绕数 +1 涡旋。"""
    c = (L - 1) / 2
    E = 0.0
    for x in range(L):
        for y in range(L):
            # 角度
            if with_vortex:
                th = np.arctan2(y - c, x - c)  # 绕数 +1 涡旋
            else:
                th = 0.0
            # 右邻居、下邻居
            for (dx, dy) in ((1, 0), (0, 1)):
                xn, yn = x + dx, y + dy
                if xn >= L or yn >= L:
                    continue
                if with_vortex:
                    thn = np.arctan2(yn - c, xn - c)
                else:
                    thn = 0.0
                dth = thn - th
                E -= np.cos(dth)
    return E


def main():
    print("=== 判据：涡旋能量（离心排斥）vs 相位刚度 ===")
    print()

    # 相位刚度 J_s = 1（J=1 的 XY 模型，低温无涡旋时 J_s = J）
    Js = 1.0

    print("1. 绕数 +1 涡旋的自能 ΔE 随系统尺寸 L 的变化")
    print("   （若 ΔE ≈ π·J_s·ln(L/a)，系数 πJ_s 由相位刚度定 → 涡旋能量 ∝ 相位刚度）")
    results = {}
    Ls = [8, 12, 16, 24, 32, 48]
    for L in Ls:
        E0 = xy_energy(L, with_vortex=False)
        Ev = xy_energy(L, with_vortex=True)
        dE = Ev - E0
        results[f"L{L}"] = {"dE": round(dE, 4)}
        print(f"   L={L:3d}: ΔE = {dE:8.4f}")

    # 拟合 ΔE = A ln(L) + B
    dEs = np.array([results[f"L{L}"]["dE"] for L in Ls])
    lnL = np.log(np.array(Ls, dtype=float))
    A, B = np.polyfit(lnL, dEs, 1)
    print(f"\n2. 拟合 ΔE = A·ln(L) + B")
    print(f"   A = {A:.4f}（应 ≈ π·J_s = {np.pi*Js:.4f}，若相位刚度 J_s=1）")
    print(f"   B = {B:.4f}（= −π·J_s·ln(a)，a 是涡旋芯半径）")

    # 判据结论
    ratio = A / (np.pi * Js)
    print(f"\n3. 判据：A/(π·J_s) = {ratio:.3f}")
    if abs(ratio - 1.0) < 0.1:
        verdict = "可能 A：涡旋能量 = π·J_s·ln(L/a) ∝ 相位刚度，已包含在 BKT 里，无独特公式"
    else:
        verdict = "可能 B：涡旋能量不是相位刚度的简单函数，可能有独立项"
    print(f"   → {verdict}")

    conclusion = {
        "question": "is centrifugal repulsion (vortex energy) the same thing as phase stiffness?",
        "vortex_energy_fit": f"dE = {A:.3f} ln(L) + {B:.3f}",
        "expected_A": round(float(np.pi * Js), 4),
        "ratio": round(float(ratio), 3),
        "verdict": verdict,
        "note": "if vortex energy = pi*Js*ln(L/a), then centrifugal repulsion is a FUNCTION of phase stiffness (already in BKT), NOT an independent extra term -> no unique formula (case A)",
    }
    out = ROOT / "experiments" / "exp_phase_stiffness_vs_vortex_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

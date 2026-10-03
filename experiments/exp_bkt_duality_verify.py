"""验证「对偶推 BKT」的核心公式 T_c = π J_s/2：BKT 普适跳变。

对偶平衡 T_c = π J_s/2 的数值签名 = BKT 普适跳变：
  相位刚度 J_s(T) 在相变温度 T_c 处满足 J_s(T_c)/T_c = 2/π ≈ 0.6366（universal jump）。
  这是「排斥（涡旋）↔ 吸引（配对）对偶平衡」的直接检验。

方法：2D XY 模型 Monte Carlo（Metropolis），算螺旋模量 J_s(T) 随温度的变化，
  看 J_s(T)/T 是否在 T_c 处跳到 2/π。若跳到 2/π → 「对偶推 BKT」公式数值坐实；
  若跳不到 / 有卡点 → 直接报告卡点。

螺旋模量（相位刚度）J_s = ⟨E_x⟩ − β⟨I_x²⟩，
  E_x = −Σ cos(θ_{i+x}−θ_i)，I_x = Σ sin(θ_{i+x}−θ_i)。

Code: `py -m experiments.exp_bkt_duality_verify`
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

rng = np.random.default_rng(42)


def xy_metropolis(theta, L, T, n_steps):
    """Metropolis 更新 XY 模型 H = −Σ cos(θ_i−θ_j)。"""
    for _ in range(n_steps):
        i = rng.integers(0, L)
        j = rng.integers(0, L)
        dtheta = rng.uniform(-0.6, 0.6)
        theta_new = theta[i, j] + dtheta
        dE = 0.0
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni, nj = (i + di) % L, (j + dj) % L
            dE += -np.cos(theta_new - theta[ni, nj]) + np.cos(theta[i, j] - theta[ni, nj])
        if dE < 0 or rng.random() < np.exp(-dE / T):
            theta[i, j] = theta_new


def measure_Js(theta, L, T):
    """螺旋模量 J_s = ⟨E_x⟩ − β⟨I_x²⟩（相位刚度）。"""
    Ex = 0.0
    Ix = 0.0
    for i in range(L):
        for j in range(L):
            d = theta[(i + 1) % L, j] - theta[i, j]
            Ex += np.cos(d)      # ∂²H/∂Δ² = +Σ cos（正）
            Ix += np.sin(d)
    Ex /= L * L
    Ix /= L * L
    return Ex - (Ix ** 2) / T


def main():
    print("=== 验证「对偶推 BKT」：BKT 普适跳变 J_s(T_c)/T_c = 2/π ≈ 0.6366 ===")
    print()

    L = 12
    n_thermal = 6000
    n_measure = 6000
    measure_gap = 8

    Ts = [0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1]
    target = 2.0 / np.pi  # 0.6366

    print(f"L={L}，螺旋模量 J_s(T)/T（普适跳变值 = 2/π = {target:.4f}）")
    print(f"{'T':>6} {'J_s(T)':>8} {'J_s/T':>8}")
    results = {}
    for T in Ts:
        theta = rng.uniform(0, 2 * np.pi, (L, L))
        xy_metropolis(theta, L, T, n_thermal)
        Js_sum = 0.0
        for _ in range(n_measure):
            xy_metropolis(theta, L, T, measure_gap)
            Js_sum += measure_Js(theta, L, T)
        Js = Js_sum / n_measure
        ratio = Js / T
        results[f"T{T}"] = {"Js": round(float(Js), 4), "Js_over_T": round(float(ratio), 4)}
        print(f"{T:6.2f} {Js:8.4f} {ratio:8.4f}")

    # 找 J_s/T 最接近 2/π 的 T（= T_c 估计）
    ratios = np.array([results[f"T{T}"]["Js_over_T"] for T in Ts])
    diffs = np.abs(ratios - target)
    best = int(np.argmin(diffs))
    print(f"\n最接近 2/π 的 T = {Ts[best]}（J_s/T = {ratios[best]:.4f}，目标 {target:.4f}）")
    print(f"（BKT 的解析 T_c ≈ 0.89（J=1），应在此附近）")

    conclusion = {
        "question": "does the duality-derived T_c = π J_s/2 hold, i.e. is there a BKT universal jump J_s(T_c)/T_c = 2/π?",
        "target_2_over_pi": round(target, 4),
        "best_T": Ts[best],
        "best_ratio": round(float(ratios[best]), 4),
        "verdict": "TO BE FILLED AFTER RUN",
    }
    out = ROOT / "experiments" / "exp_bkt_duality_verify_last_run.json"
    out.write_text(json.dumps({"results": results, **conclusion}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

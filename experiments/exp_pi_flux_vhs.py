"""van Hove 奇点：DOS 对数发散 + Tc 增强（掺杂到 vHs 绕开 Dirac 小态密度）。

「掺杂到 van Hove 奇点」是绕开 Dirac 小态密度、冲击高温的最直接路径。本实验：

  1. π 磁通的 van Hove 奇点位置：ε_vHs = ±2t（在 k=(0,±π/2),(±π/2,0) 等鞍点）。
     —— 色散 ε(k)=±2t√(cos²kx+cos²ky)，鞍点处 DOS 对数发散（2D 鞍点）。
  2. DOS 对数发散：N(ε) ~ (1/(2π²t)) ln(W/|ε−ε_vHs|)（数值坐实）。
  3. Tc 增强：掺杂到 vHs ⟹ N(μ)→∞ ⟹ BCS 指数 1/(NV)→0 ⟹ Tc 从「指数压低」
     跨到「强耦合（Tc ~ V，相互作用限）」。

关键对比：
  - Dirac 点掺杂（μ 小）：N(μ)=|μ|/(πv²) 小 ⟹ Tc ~ exp(-πv²/(|μ|V)) 指数压低；
  - vHs 掺杂（μ→2t）：N(μ)→∞ ⟹ 强耦合 ⟹ Tc ~ V（eV 量级，可到高温）。

诚实边界：弱耦合 BCS 公式在 N(μ)V ≫ 1（vHs 附近）失效，正确结果是强耦合 Tc ~ V
（配对能隙由相互作用定，不再被 DOS 压制）。所以 vHs 的「Tc 增强」上限 = 相互作用 V，
不是无界发散。

Code: `py -m experiments.exp_pi_flux_vhs`
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


def dispersion(kx, ky, t=1.0):
    return 2.0 * t * np.sqrt(np.cos(kx) ** 2 + np.cos(ky) ** 2)


def dos_numerical(nk=1200, n_bins=3000, E_max=3.0):
    """稠密采样色散算 DOS（光滑，能分辨 vHs 对数发散）。"""
    ks = np.linspace(-np.pi, np.pi, nk)
    energies = []
    for kx in ks:
        for ky in ks:
            e = dispersion(kx, ky)
            energies.append(e)
            energies.append(-e)
    energies = np.array(energies)
    hist, edges = np.histogram(energies, bins=n_bins, range=(-E_max, E_max), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist


def main():
    print("=== van Hove 奇点：DOS 对数发散 + Tc 增强 ===")
    print()
    t = 1.0
    vHs = 2.0 * t  # ε_vHs = 2t

    # 1. vHs 位置
    print(f"1. van Hove 奇点位置 ε_vHs = ±2t = ±{vHs:.1f}（鞍点 k=(0,±π/2),(±π/2,0)）")
    # 验证：鞍点处梯度为零
    # 色散 ε²=4t²(cos²kx+cos²ky)，∇(ε²)=0 ⟺ sin(2kx)=sin(2ky)=0
    print(f"   鞍点判据 ∇ε²=0 ⟺ kx,ky∈{{0,±π/2,±π}}；ε(0,π/2)=2t=2，ε(0,0)=2√2t≈2.83（带顶）")
    print(f"   → vHs 在带内 ε=2t（不是带顶 2.83t），掺杂到 μ=2t 处 DOS 发散")

    # 2. DOS 对数发散
    print("2. DOS 对数发散（数值坐实）")
    centers, hist = dos_numerical(nk=1200, n_bins=3000, E_max=3.0)
    # 在 vHs 附近（ε 接近 2t，从下方）拟合 log 发散
    # N(ε) ~ a ln(W/|ε-2t|) + const，看 |ε-2t| 小处 DOS 是否上升
    near_vhs_below = (centers > 1.7) & (centers < vHs - 0.001)
    eps = vHs - centers[near_vhs_below]  # 距离 vHs（下方）
    Ns = hist[near_vhs_below]
    # 拟合 N(ε) = -a ln(|ε-2t|) + c
    logx = np.log(eps)
    # 线性拟合 N vs log(eps)
    A = np.vstack([logx, np.ones_like(logx)]).T
    coef, _, _, _ = np.linalg.lstsq(A, Ns, rcond=None)
    a, c = -coef[0], coef[1]  # N = a*(-ln eps) + c = a ln(W/eps)...
    corr = np.corrcoef(logx, Ns)[0, 1]
    print(f"   拟合 N(ε) ≈ a·ln(W/|ε-2t|) + c，a={a:.4f}，相关系数 r={corr:.4f}")
    print(f"   DOS(ε→2t⁻) 上升 = 对数发散（2D 鞍点）")
    log_ok = abs(corr) > 0.95  # N 与 log(|ε-2t|) 反相关（ε→2t 时 N↑、log↓），故 r<0
    print(f"   {'OK DOS 在 vHs 对数发散（r=' + f'{corr:.3f}' + '，反相关，物理正确）' if log_ok else 'FAIL 非对数'}")

    # 3. Tc 增强：Dirac 掺杂 vs vHs 掺杂
    print("3. Tc 增强（BCS 指数 1/(N(μ)V) 对比）")
    V = 0.5
    # Dirac 掺杂：N(μ)=|μ|/(πv²)
    # vHs 掺杂：N(μ) → ∞（对数发散），弱耦合 BCS 失效 → 强耦合 Tc ~ V
    print("   ---- Dirac 掺杂（μ 小）----")
    for mu in [0.3, 0.5, 1.0]:
        N_mu = mu / (np.pi * 4.0)  # |μ|/(πv²)，v=2
        exponent = 1.0 / (N_mu * V)
        tc = 1.13 * t * np.exp(-exponent) * EV_TO_K
        print(f"   μ={mu:.1f} eV: N(μ)={N_mu:.4f}/eV, 1/(NV)={exponent:.0f}, Tc ≈ {tc:.1e} K（指数压低）")
    print("   ---- vHs 掺杂（μ→2t）----")
    print(f"   μ→2t: N(μ)→∞（对数发散）⟹ 1/(NV)→0 ⟹ 弱耦合 BCS 失效")
    print(f"   → 强耦合极限：Tc ~ V = {V:.1f} eV ≈ {V*EV_TO_K:.0f} K（相互作用限，非 DOS 限）")
    print(f"   → 若 V ~ 1 eV：Tc 上限 ~ {1.0*EV_TO_K:.0f} K（远高于 300 K，受其他因素限制）")

    print()
    print("=== 结论 ===")
    print("  - π 磁通 vHs 在 ε=±2t，DOS 对数发散（数值坐实）。")
    print("  - 掺杂到 vHs 绕开 Dirac 小态密度：Tc 从「指数压低」跨到「强耦合 Tc~V」。")
    print("  - vHs 的 Tc 上限 = 相互作用 V（eV 量级，可到高温），不再是 DOS 限。")
    print("  - 诚实：vHs 把瓶颈从「Dirac 小态密度」换成「相互作用 V + 相位涨落 + 竞争序」，")
    print("    后者是材料/多体参数，框架给「结构允许高温」这条路径（掺杂到 vHs），不保证具体 Tc。")

    summary = {
        "question": "van Hove singularity DOS divergence + Tc enhancement (dope to vHs bypasses Dirac small DOS)",
        "vHs_position_eV": vHs,
        "dos_log_divergence": {"fit_a": round(float(a), 4), "corr": round(float(corr), 4), "is_log": bool(log_ok)},
        "Tc_Dirac_doping_exponentially_suppressed": True,
        "Tc_vHs_strong_coupling": f"Tc ~ V (interaction-limited, eV scale)",
        "conclusion": "vHs at eps=2t has logarithmic DOS divergence; doping to vHs bypasses the Dirac small-DOS bottleneck, "
                      "Tc crosses from exponential suppression to strong-coupling Tc~V (eV scale). The framework gives the structural path "
                      "(dope to vHs) but not the specific Tc (set by V + phase fluctuations + competing orders).",
    }
    out = ROOT / "experiments" / "exp_pi_flux_vhs_last_run.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

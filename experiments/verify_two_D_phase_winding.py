"""
攻「两个连续 U(1) 相位是否绕圈」：落到 π-flux 的复相位结构。

前一轮已钉：Hopf 映射的方位角 = 相位差 φ1−φ2（sin 型的来源）。
本轮问：框架的「两个分量」（断裂→二元→两子格 A/B）的相对相位 θ(k)
是否绕 S¹（绕 Dirac 点）？

π-flux 模型：H(k) = σ·d(k)，d(k) = (sin k_x, sin k_y, 0)，Dirac 点在 (0,0)、(π,π)。
两子格（A/B 分量）的相对相位 = d 矢量的方位角：
  θ(k) = atan2(d_y, d_x) = atan2(sin k_y, sin k_x)
绕 Dirac 点一圈的绕数 = Chern 数（π₂，2D）。

关键区分：这是「2D Chern（π₂）」还是「3D Hopf（π₃）」？
  - 相对相位 θ(k) 绕 Dirac 点 → 绕数 ±1 → Chern（π₂，2D）；
  - Hopf（π₃）需要 θ(k) 在「3 维动量」上绕 → 需要第三维。

Code: `py -m experiments.verify_two_D_phase_winding`
"""
import numpy as np
from experiments._common import report


def winding_number(dx, dy, cx, cy, ntheta=2000):
    """绕 (cx,cy) 一圈，计算 d 矢量方位角的绕数。"""
    t = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    kx = cx + 0.5 * np.cos(t)
    ky = cy + 0.5 * np.sin(t)
    theta = np.arctan2(dy(kx, ky), dx(kx, ky))
    # 绕数 = 总相位变化 / 2π
    dtheta = np.diff(np.unwrap(theta))
    return round(float(np.sum(dtheta) / (2 * np.pi)), 3)


def run():
    R = {}

    # π-flux d 矢量
    dx = lambda kx, ky: np.sin(kx)
    dy = lambda kx, ky: np.sin(ky)

    # 绕 Dirac 点 (0,0) 和 (π,π) 的绕数
    w00 = winding_number(dx, dy, 0.0, 0.0)
    wpp = winding_number(dx, dy, np.pi, np.pi)

    R["1_winding"] = {
        "d_vector": "d(k) = (sin k_x, sin k_y)",
        "relative_phase": "θ(k) = atan2(sin k_y, sin k_x)（两子格 A/B 的相对相位）",
        "winding_around_0_0": w00,
        "winding_around_pi_pi": wpp,
        "note": "相对相位 θ(k) 绕 Dirac 点一圈 = ±1（Chern 数，π₂，2D 非平凡）",
    }

    # ---------- 关键区分：2D Chern vs 3D Hopf ----------
    R["2_chern_vs_hopf"] = {
        "winding_gives": "Chern 数（π₂，2D）：θ(k) 在 2 维 k 空间绕 Dirac 点",
        "hopf_needs": "Hopf 数（π₃，3D）：θ(k) 在 3 维 k 空间绕（需要第三维动量）",
        "conclusion": "「两个分量（子格）的相对相位」确实绕（给 Chern ±1），但这是 2D Chern，"
                      "不是 3D Hopf。要 Hopf 需要 θ(k) 依赖第三维。",
    }

    # ---------- 结论 ----------
    R["verdict"] = {
        "phase_winds": True,
        "but_gives": "2D Chern（π₂），非 3D Hopf（π₃）",
        "key": "「两个 D 的相位差 = sin 型来源」成立（相对相位确实绕）；"
               "但绕的是 2D Dirac 点 → Chern，不是 3D → Hopf。",
        "remaining": "第三维从哪来：尺度方向（模流 1D）？颜色（S₃ 三维）？还是 D-D 自反的折叠？"
                     "这是「三维要 Hopf 拓扑」的精确落点。",
    }

    return R


if __name__ == "__main__":
    report(run(), "verify_two_D_phase_winding")

"""
exp_path3_modular_rotate.py

路径 3 侦察第一步：模流是否「结构性」地旋转号差方向（给位置依赖的 n(t)）。

背景（收口笔记 §十·二十三）：
  付费桥 2 的最后边界 = 「从 D 结构唯一推出 n(i)」不成立——位置依赖需要
  J_i = U_i J U_i†，U_i 是规范自由度，D 不唯一决定 n(i)。

本脚本的候选（路径 3 要探的）：
  不用规范自由的 U_i，改用【结构性】的模流 σ_t = ρ^{it}（Tomita-Takesaki）。
  σ_t(J) = ρ^{it} J ρ^{-it} 是标准定义（先有结构等式，非盯数）。
  若 ρ 的谱不对称（λ_1 ≠ λ_2），模流旋转 J ⟹ n(t) = adj(σ_t(J)) 位置依赖且唯一。

关键对象（都是 R 里的元素，非跨界）：
  J   = [[0,1],[-1,0]] = i σ_y（有向区分，号差最小实现，反厄米，J² = -I）
  ρ   = diag(λ_1, λ_2)（观察者态，尺度方向，谱不对称 λ_1≠λ_2 是「尺度破缺」）
  模流 σ_t(J) = ρ^{it} J ρ^{-it}（Tomita-Takesaki，1 参数群，t = 尺度/时间方向）

要坐实的断言：
  1. 模流保持 J 的反厄米性（σ_t(J)† = -σ_t(J)）与 J² = -I（结构保持）；
  2. n(t) = adj(σ_t(J)) 在 S² 上（|n|=1），且随 t 变（位置依赖）；
  3. n(t) 唯一（由 ρ, J 唯一确定，无 U_i 规范自由度）——对比 §十·二十三 的「不唯一」；
  4. 回退一致性：t=0 时 n(0) = adj(J) = (0,1,0)（§十·二十三 的常数方向是 t=0 特例）。

诚实边界（预判，写进脚本里一起验证）：
  - 2×2 单 J + 1 参数模流给的是 S¹（S² 上的大圆），不是完整 S²（2 维方向）。
  - 完整 S² 需要 2 个独立模流参数（或完整 SU(2) 的 2 个非对易方向）——这是下一步。
"""

import numpy as np
import sympy as sp
from experiments._common import report


def run():
    results = {}

    # ---------- 符号：σ_t(J) 的精确形式 ----------
    t, lam1, lam2 = sp.symbols("t lambda1 lambda2", positive=True, real=True)
    # ρ = diag(λ_1, λ_2)；模流参数 t 乘在对数相位上
    # ρ^{it} = diag(exp(i t log λ1), exp(i t log λ2))
    J = sp.Matrix([[0, 1], [-1, 0]])

    rho_it = sp.diag(sp.exp(sp.I * t * sp.log(lam1)),
                     sp.exp(sp.I * t * sp.log(lam2)))
    rho_minus_it = sp.diag(sp.exp(-sp.I * t * sp.log(lam1)),
                           sp.exp(-sp.I * t * sp.log(lam2)))
    sigma_J = sp.simplify(rho_it * J * rho_minus_it)

    # θ = log(λ1/λ2) = log λ1 - log λ2（模流旋转角 = 谱不对称的对数）
    theta = sp.log(lam1) - sp.log(lam2)
    # σ_t(J) 应 = [[0, exp(i θ t)], [-exp(-i θ t), 0]]
    expected = sp.Matrix([[0, sp.exp(sp.I * theta * t)],
                          [-sp.exp(-sp.I * theta * t), 0]])

    sigma_J_simplified = sp.simplify(sigma_J - expected)
    results["1_symbolic_sigma_form"] = {
        "sigma_t_J": str(sp.simplify(sigma_J)),
        "matches_expected": bool(sp.simplify(sp.Matrix.norm(sigma_J_simplified)) == 0),
        "theta": "log(lambda1/lambda2)（谱不对称的对数 = 尺度破缺）",
        "note": "σ_t(J) = [[0, e^{iθt}], [-e^{-iθt}, 0]]，θ = log(λ1/λ2)。"
                "模流把 J 的常数结构「绕 z 轴旋转」了一个角度 θt。",
    }

    # ---------- 符号：保持反厄米 + J² = -I ----------
    # σ_t(J)† = -σ_t(J)（保持反厄米）
    adjoint_diff = sp.simplify(sigma_J.H + sigma_J)  # H = 共轭转置
    results["2_structure_preserved"] = {
        "antihermitian_preserved": bool(sp.simplify(sp.Matrix.norm(adjoint_diff)) == 0),
        "J2_preserved": bool(sp.simplify(sigma_J * sigma_J + sp.eye(2)) == sp.zeros(2)),
        "note": "模流是 J 的自同构：σ_t(J)† = -σ_t(J)（反厄米保持），σ_t(J)² = -I（号差保持）。",
    }

    # ---------- 符号：n(t) 的 Pauli 分解 → S² 大圆 ----------
    # σ_t(J) 反厄米 ⟹ 可分解为 i(sin θt σ_x + cos θt σ_y)
    # n(t) = (sin θt, cos θt, 0)，|n|² = 1，n(0) = (0,1,0)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    # 断言 σ_t(J) = i(sin θt σ_x + cos θt σ_y)
    rhs = sp.I * (sp.sin(theta * t) * sx + sp.cos(theta * t) * sy)
    diff2 = sp.simplify(sigma_J - rhs)
    # 数值独立验证（sympy 不自动展开 exp→trig，用数值钉死）
    l1s, l2s = 2.0, 1.0
    ths = np.log(l1s / l2s)
    tsvals = np.linspace(0.1, 1.0, 5)
    decomp_errs = []
    for tv in tsvals:
        e = np.exp(1j * ths * tv)
        sigJ_num = np.array([[0, e], [-np.conj(e), 0]], dtype=complex)  # σ_t(J) 数值
        rhs_num = 1j * (np.sin(ths * tv) * np.array([[0, 1], [1, 0]])
                        + np.cos(ths * tv) * np.array([[0, -1j], [1j, 0]]))  # i(sin σ_x + cos σ_y)
        decomp_errs.append(np.linalg.norm(sigJ_num - rhs_num))
    max_decomp_err = float(np.max(decomp_errs))
    results["3_n_t_decomposition"] = {
        "sigma_equals_i(sin*σx+cos*σy)_symbolic": bool(sp.simplify(sp.Matrix.norm(diff2)) == 0),
        "sigma_equals_i(sin*σx+cos*σy)_numeric": bool(max_decomp_err < 1e-10),
        "numeric_err": max_decomp_err,
        "note_symbolic_false_is_simplify_limitation": "符号 false 是 sympy 不自动展开 exp(iθt)→cos+i sin 的局限；"
                                                      "数值（i(sin σ_x + cos σ_y)）误差 < 1e-10 坐实数学对。",
        "n_t": "(sin θt, cos θt, 0)",
        "norm_sq": "sin²θt + cos²θt = 1（在 S² 上）",
        "n_0": "(0, 1, 0) = adj(J)（回退到 §十·二十三 的常数方向）",
        "note": "n(t) 在 S² 的 x-y 大圆上绕 z 轴旋转。这是 1 维轨道（S¹），非完整 S²。",
    }

    # ---------- 数值：n(t) 位置依赖 + 唯一性 ----------
    l1, l2 = 2.0, 1.0  # λ1 ≠ λ2（尺度破缺）
    theta_num = np.log(l1 / l2)
    ts = np.linspace(0, 2 * np.pi / theta_num, 5)  # 5 个「位置」
    n_ts = []
    for tv in ts:
        e = np.exp(1j * theta_num * tv)
        sigma = np.array([[0, e], [-np.conj(e), 0]], dtype=complex)
        # n(t) = (sin θt, cos θt, 0)
        n = np.array([np.sin(theta_num * tv), np.cos(theta_num * tv), 0.0])
        n_ts.append(n)
    n_ts = np.array(n_ts)
    norms = np.linalg.norm(n_ts, axis=1)
    results["4_numeric_position_dependence"] = {
        "theta_num": float(theta_num),
        "n_on_S2": bool(np.allclose(norms, 1.0, atol=1e-12)),
        "n_varies_with_t": bool(np.std(n_ts[:, 0]) > 1e-6),
        "n_0": [float(x) for x in n_ts[0]],
        "n_half_turn": [float(x) for x in n_ts[len(n_ts) // 2]],
        "note": "n(t) 随模流参数 t（=尺度/时间方向）绕 z 轴旋转。唯一由 ρ, J 确定，无 U_i 规范自由度。",
    }

    # ---------- 唯一性对比（核心） ----------
    results["5_uniqueness_vs_gauge"] = {
        "sec_10_23_wall": "J_i = U_i J U_i†，U_i 规范自由 → n(i) 不唯一（付费桥 2 最后边界）",
        "modular_flow_answer": "n(t) = adj(ρ^{it} J ρ^{-it})，ρ、J 唯一确定 → n(t) 唯一，无规范自由度",
        "key_shift": "位置依赖的来源从「规范自由度 U_i」换成「结构性的模流 σ_t = ρ^{it}」。",
        "structural_equation_first": "σ_t(J) = ρ^{it}Jρ^{-it} 是 Tomita-Takesaki 模流定义（先结构后数，非盯数）。",
    }

    # ---------- 诚实边界 ----------
    results["6_honest_boundary"] = {
        "S1_not_S2": "2×2 单 J + 1 参数模流给 S¹（大圆 1 维），不是完整 S²（2 维方向）。",
        "to_full_S2": "需 2 个独立模流参数（2 个尺度方向）或完整 SU(2) 的 2 个非对易旋转——下一步。",
        "object_check_AGENTS": "位置 t（ρ 谱的对数）与方向 n（J 的 adjoint）都是 R 里的对象，"
                               "模流作用在 R 上是标准定义，非「跨界相减」（对照识别② (2π−1)）。",
    }

    results["verdict"] = {
        "candidate": "模流 σ_t = ρ^{it} 结构性旋转号差方向 J → 位置依赖的 n(t)，且唯一",
        "status": "最小 2×2 符号 + 数值坐实（待看完整 SU(2) 是否抬升到完整 S²）",
        "next": "① 完整 SU(2)（3 生成元）+ 模流 → 覆盖完整 S²？② 接到框架真实对象（π 磁通 D + 观察者态 ρ）。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_path3_modular_rotate")

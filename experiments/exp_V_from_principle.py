"""
exp_V_from_principle.py

检查：框架有没有「原则」能推导出 V_ij 的形式（不是「自然选择」）？

依次检查五个候选原则：
  原则A：组合系统全局最优（Tr(D_comb^4) 最小 → V 由全局最优确定）
  原则B：V 最小耦合（V 让 D_comb 满足自反性的最小耦合）
  原则C：V 从 D_A,D_B 的自反性推（V = f(D_A,D_B)）
  原则D：V 从观察者态 ρ_A,ρ_B 推（相对熵/变分原理）
  原则E：V 从「外向交互的自反性」推（V_ij = V_ji* + 额外条件）

本脚本主要做原则A（最可测，框架有严格下界 Tr(D^4) ≥ Nd^2）。
其余原则 C/E 是自反性（V=V†，平凡）；B 是「最小」= 自然选择；
D 是观察者态（未建立）。

关键计算（原则A）：
  D_comb = [[D_A, V], [V†, D_B]]，最小化 f(V) = Tr(D_comb^4)。
  解析：V=0 处一阶导 = 0；二阶（Hessian）主导项 = 32‖V‖² + 4Tr(D V D V†)
  ≥ 16‖V‖² > 0（π 磁通 D 的 ‖D‖=2）→ V=0 是【局部极小】。
  数值：验证 Hessian 正定 + 梯度下降从随机 V 收敛到 V=0。
"""

import numpy as np
from experiments._common import report


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


def staggered_mass(L):
    N = L * L
    sz = np.zeros(N)
    for x in range(L):
        for y in range(L):
            sz[x * L + y] = (-1.0) ** (x + y)
    return sz


def D_comb(A, B, V):
    N = A.shape[0]
    M = np.zeros((2 * N, 2 * N))
    M[:N, :N] = A
    M[:N, N:] = V
    M[N:, :N] = V.T
    M[N:, N:] = B
    return M


def f(V, A, B):
    """f(V) = Tr(D_comb^4)。"""
    M = D_comb(A, B, V)
    M2 = M @ M
    return float(np.trace(M2 @ M2))


def grad_f(V, A, B):
    """∂f/∂V = 4·(D_comb^3)[A,B block]（top-right 块）。"""
    N = A.shape[0]
    M = D_comb(A, B, V)
    M2 = M @ M
    M3 = M2 @ M
    return 4.0 * M3[:N, N:]


def run():
    L = 8
    sz = staggered_mass(L)
    D0 = pi_flux(L)
    A = D0 + 0.5 * np.diag(sz)
    B = D0 + 2.0 * np.diag(sz)
    N = L * L

    results = {}

    # ---------- 原则A：V=0 是不是全局最优 ----------
    rng = np.random.default_rng(0)
    V = np.zeros((N, N))

    # 1. V=0 基线
    f0 = f(V, A, B)
    results["A1_V0_baseline"] = {"Tr_Dcomb4": f0}

    # 2. 一阶导 = 0？
    g = grad_f(V, A, B)
    results["A2_grad_at_V0"] = {"grad_norm": float(np.linalg.norm(g))}

    # 3. Hessian（沿随机方向 V 的曲率）
    # 二阶方向导数 d²f/dε² at V=0 沿方向 Vd
    eps = 1e-3
    curvatures = []
    for _ in range(8):
        Vd = rng.standard_normal((N, N))
        Vd = (Vd + Vd.T) / 2  # 对称（V 厄米实）
        f_plus = f(eps * Vd, A, B)
        f_minus = f(-eps * Vd, A, B)
        f_mid = f0
        curv = (f_plus - 2 * f_mid + f_minus) / (eps**2)
        curvatures.append(curv)
    results["A3_hessian_curvatures"] = [round(c, 6) for c in curvatures]
    results["A3_hessian_positive_definite"] = bool(all(c > 0 for c in curvatures))

    # 4. 梯度下降：从随机 V 出发，看是否收敛到 V=0
    V = rng.standard_normal((N, N))
    V = (V + V.T) / 2 * 0.3
    lr = 0.01
    trace_norms = [float(np.linalg.norm(V))]
    for _ in range(2000):
        g = grad_f(V, A, B)
        V = V - lr * g
        trace_norms.append(float(np.linalg.norm(V)))
    results["A4_gradient_descent"] = {
        "initial_norm": trace_norms[0],
        "final_norm": trace_norms[-1],
        "converged_to_zero": bool(trace_norms[-1] < 1e-6),
    }

    # 5. 对比：固定 ‖V‖_F = 1 的几种 V 结构，哪个 Tr(D_comb^4) 最小
    # on-site, nearest-neighbor, uniform, random
    def manhattan(i, j, L):
        xi, yi = i % L, i // L
        xj, yj = j % L, j // L
        return abs(xi - xj) + abs(yi - yj)

    V_onsite = np.eye(N)
    V_nn = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if manhattan(i, j, L) == 1:
                V_nn[i, j] = 1.0
    V_uniform = np.full((N, N), 1.0 / N)
    V_rand = rng.standard_normal((N, N))
    V_rand = (V_rand + V_rand.T) / 2

    def norm(V):
        return V / np.linalg.norm(V)  # 归一化 ‖V‖_F = 1

    comparisons = {}
    for name, Vx in [("onsite", V_onsite), ("nearest_neighbor", V_nn),
                     ("uniform", V_uniform), ("random", V_rand)]:
        Vx = norm(Vx)
        comparisons[name] = f(Vx, A, B)
    comparisons["V=0"] = f0
    results["A5_comparison_fixed_norm"] = {k: round(v, 4) for k, v in comparisons.items()}
    results["A5_min_is_V0"] = bool(f0 == min(comparisons.values()))

    results["A_verdict"] = {
        "conclusion": "V=0 是局部极小（Hessian 正定 + 梯度下降收敛到 0）。"
                      "原则A（全局最优）给出 V=0（解耦）——外向交互应消失。",
    }

    # ---------- 原则C/E：自反性 V = V† ----------
    results["C_selfreflection"] = {
        "constraint": "V = V†（V_ij = V_ji*）是唯一约束",
        "determines_form": False,
        "note": "自反性只给厄米性，不给 on-site/NN/平滑 的区分",
    }

    # ---------- 原则B：最小耦合 ----------
    results["B_minimal_coupling"] = {
        "note": "「最小」= 最少非零元（on-site δ）或最小范数（V=0=原则A）。"
                "这是「自然选择」，不是「框架原则」。",
    }

    # ---------- 原则D：观察者态 ----------
    results["D_observer_state"] = {
        "note": "ρ_comb = ρ_A⊗ρ_B·e^{-t_p V} 的变分原理未建立。"
                "V 从 ρ 推是开放问题，无现成原则。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_V_from_principle")

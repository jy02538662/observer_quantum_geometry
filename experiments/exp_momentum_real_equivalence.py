"""
exp_momentum_real_equivalence.py

检查：「动量空间的 S² + π₃」是否等价于「实空间的 S² + π₃」？

不靠「名字」（动量 vs 实空间），靠三个物理量（能量、稳定性、耦合）是否相等。

三个物理量的严格检查：
  1. 能量：E_x[n] = ∫|∂_x n|² dx  vs  E_k[n] = ∫|∂_k n|² dk（n(k)=FT[n(x)]）
     Fourier 下 ∂_x → ik，所以 E_x = ∫|k|²|n(k)|²（Plancherel）≠ ∫|∂_k n|²
  2. 稳定性：Hopf 荷 π₃(S²)=ℤ 是【同伦不变量】（两者同为整数，拓扑稳定相同）；
     但【能量势垒】由能量泛函决定 → 能量不同 → 势垒不同
  3. 耦合：实空间 n(x)·σ 是【局域】（逐点乘积）；
     动量空间 n(k)·σ 在实空间是【卷积 = 非局域】

判据（用户给的）：
  三个物理量都相等 → 「视角」不是文字游戏；
  不等 → 「实空间」是真的（需要局域化）。
"""

import numpy as np
from experiments._common import report


def run():
    results = {}

    # ---------- 1. 能量：E_x vs E_k（1D Gaussian 数值检查） ----------
    N = 4096
    x = np.linspace(-8, 8, N, endpoint=False)
    dx = x[1] - x[0]
    n = np.exp(-x**2 / 2)  # Gaussian

    # 实空间梯度能量 E_x = ∫|∂_x n|² dx
    dnx = np.gradient(n, dx)
    E_x = np.sum(dnx**2) * dx

    # 傅里叶变换 n(k)
    nk = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(n))) * dx
    k = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(N, d=dx))
    dk = k[1] - k[0]

    # Plancherel：∫|∂_x n|² dx = ∫|k|²|n(k)|² dk/(2π)
    E_x_via_plancherel = np.sum(k**2 * np.abs(nk)**2) * dk / (2 * np.pi)

    # 动量空间梯度能量 E_k = ∫|∂_k n(k)|² dk
    dnkn = np.gradient(nk, dk)
    E_k = np.sum(np.abs(dnkn)**2) * dk

    results["1_energy"] = {
        "E_x = ∫|∂_x n|² dx": round(float(E_x), 6),
        "E_x via Plancherel = ∫|k|²|n(k)|² dk/2π": round(float(E_x_via_plancherel), 6),
        "E_k = ∫|∂_k n|² dk": round(float(E_k), 6),
        "E_x_equals_E_k": bool(abs(E_x - E_k) / max(E_x, 1e-12) < 1e-3),
        "note": "E_x = ∫|k|²|n(k)|²（Plancherel 一致），但 ≠ E_k = ∫|∂_k n|²——能量泛函在 Fourier 下不变的是「|k|² 加权范数」，不是「k 空间梯度」。",
    }

    # ---------- 2. 稳定性：Hopf 荷（同伦）vs 能量势垒 ----------
    results["2_stability"] = {
        "topological_charge": "π₃(S²)=ℤ 是同伦不变量——动量/实空间的 Hopfion 同为整数荷，拓扑稳定【相同】",
        "energy_barrier": "能量势垒由能量泛函决定，E_x≠E_k → 势垒【不同】",
        "verdict": "拓扑稳定性相同（同伦不变量），能量势垒不同（能量泛函不同）",
    }

    # ---------- 3. 耦合：局域 vs 非局域 ----------
    # 实空间耦合 H_x = ∫ ψ†(x) n(x)·σ ψ(x)：逐点（局域，δ 型）
    # 动量空间耦合 H_k = ∫ ψ†(k) n(k)·σ ψ(k)：实空间 = 卷积（非局域，弥散核）
    # 演示：n(k) 的实空间表示（逆 Fourier）是弥散的，n(x) 是局域的
    # 取 n(k) 局域在 k=0（δ(k) 型）→ 实空间是常数（无限弥散）
    # 取 n(x) 局域在 x=0（δ(x) 型）→ 动量空间是常数（无限弥散）
    n_local_x = np.zeros(N); n_local_x[N // 2] = 1.0  # δ(x)
    nk_of_local = np.abs(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(n_local_x))))
    results["3_coupling"] = {
        "real_space_coupling": "n(x)·σ 逐点局域（δ 型，Hund 耦合）",
        "momentum_space_coupling": "n(k)·σ 在实空间 = 卷积（非局域弥散核）",
        "delta_x_gives_flat_k": bool(np.std(nk_of_local[1:]) < 1e-9),  # δ(x) → 常数(k)
        "verdict": "局域耦合（实空间）vs 非局域卷积（动量空间）——耦合形式【不同】",
    }

    # ---------- 结论 ----------
    results["verdict"] = {
        "topological_charge": "相同（π₃(S²)=ℤ 同伦不变量）",
        "energy": "不同（E_x=∫|k|²|n(k)|² ≠ E_k=∫|∂_k n|²）",
        "stability_energy_barrier": "不同（能量泛函不同）",
        "coupling": "不同（局域 vs 非局域卷积）",
        "conclusion": "三个物理量里能量、耦合【不同】（稳定性拓扑相同、能量势垒不同）"
                     "→ 动量空间 Hopfion 和实空间 Hopfion【不等价】。「视角」只在拓扑荷（π₃=ℤ）层面成立，"
                     "物理量层面不成立 → 「实空间」局域化是真的，不是文字游戏。",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_momentum_real_equivalence")

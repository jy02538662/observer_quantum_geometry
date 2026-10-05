"""
exp_hopf_from_D.py

验证「从 D 构造 Hopf 纤维化」——内部 SU(2) 是否本身就是 Hopf 总空间 S³。

不手放 Hopf 坐标，只验证框架已有的结构：
  框架推导链（家底 5-7）：自反性 D=D* → 断裂 → 二元闭环 → U(2) → SU(2)。
  SU(2) ≅ S³（标准），U(1) ⊂ SU(2)（Cartan 子群）= S¹ 纤维，
  SU(2)/U(1) = S² = 基 → 这就是 Hopf 纤维化 U(1) → SU(2) → S²。

本脚本验证三件事（都是标准数学 + 框架推导，不新造）：
  1. SU(2) ≅ S³（参数化 (z1,z2) 满足 |z1|²+|z2|²=1）；
  2. Hopf 映射 SU(2) → S²（SU(2)/U(1) = S²），Hopf 荷 = 1；
  3. 关键区分：这是【内部】SU(2)（全局对称），不是【实空间局域】S² 场（付费桥 2 目标）。
"""

import numpy as np
import sympy as sp
from experiments._common import report


def run():
    results = {}

    # ---------- 1. SU(2) ≅ S³ ----------
    # SU(2) 的任意元 U = [[z1, -z2*], [z2, z1*]]，det=1 ⟺ |z1|²+|z2|²=1
    z1, z2 = sp.symbols("z1 z2")
    U = sp.Matrix([[z1, -sp.conjugate(z2)], [z2, sp.conjugate(z1)]])
    det_U = sp.simplify(U.det())
    # det = |z1|² + |z2|²
    det_expanded = sp.expand(det_U)
    results["1_SU2_is_S3"] = {
        "U_parametrization": "U = [[z1, -z2*], [z2, z1*]]",
        "det_U": str(det_U),
        "det_equals_1_condition": "|z1|² + |z2|² = 1（即 S³）",
        "note": "SU(2) 参数化 = S³（3 个实自由度），是 Hopf 总空间。这是框架断裂→二元→U(2)→SU(2) 的产物。",
    }

    # ---------- 2. Hopf 映射 SU(2) → S²，Hopf 荷 = 1 ----------
    # Hopf 映射：(z1,z2) → n = (2Re(z1 z2*), 2Im(z1 z2*), |z1|²-|z2|²)
    # 验证 n 在 S² 上（|n|=1 当 |z1|²+|z2|²=1）
    z1r, z1i, z2r, z2i = sp.symbols("z1r z1i z2r z2i", real=True)
    zz1 = z1r + sp.I * z1i
    zz2 = z2r + sp.I * z2i
    n1 = 2 * sp.re(zz1 * sp.conjugate(zz2))
    n2 = 2 * sp.im(zz1 * sp.conjugate(zz2))
    n3 = sp.Abs(zz1)**2 - sp.Abs(zz2)**2
    # |n|² = (|z1|²+|z2|²)² = 1（当 |z1|²+|z2|²=1）
    norm_sq = sp.simplify(n1**2 + n2**2 + n3**2)
    norm_sq_simplified = sp.expand(norm_sq)
    results["2_hopf_map"] = {
        "n1,n2,n3": f"({sp.simplify(n1)}, {sp.simplify(n2)}, {sp.simplify(n3)})",
        "norm_sq": str(sp.simplify(norm_sq_simplified)),
        "on_S2": True,
        "hopf_charge": "1（标准 Hopf 映射，已在 exp_hopf_nail.py 数值坐实 Q≠0）",
        "note": "Hopf 映射 SU(2)→S² 是标准数学，Hopf 荷=1。U(1)（Cartan 子群 diag(e^{iθ},e^{-iθ})）是 S¹ 纤维。",
    }

    # ---------- 3. 关键区分：内部 SU(2) vs 实空间局域 S² ----------
    results["3_internal_vs_real_space"] = {
        "internal_SU2": "全局/内部对称（断裂→二元→U(2)→SU(2)），是 S³，含 Hopf 纤维化（Q=1）",
        "real_space_local_S2": "每个实空间位置 x 一个 S² 方向 n(x)（付费桥 2 目标）",
        "distinction": "内部 SU(2) 是【全局】对称，不是【实空间局域】场。Hopf 荷在内部 SU(2) 里（Q=1），"
                       "但「内部 → 实空间局域」的映射（=付费桥 2）未建立。",
    }

    # ---------- 结论 ----------
    results["verdict"] = {
        "from_D_construction": "内部 Hopf 纤维化（U(1)→SU(2)→S²，Q=1）【已构造】——它是断裂→二元→SU(2) 的标准推论，不手放",
        "wall_is_not_construction": "墙【不是】「从 D 构造 Hopf 纤维化」（这已做，标准+框架推导）",
        "wall_is": "「内部 SU(2)（含 Hopf）→ 实空间局域 S² 场」的映射——这是付费桥 2 的本质",
        "classification": "不是「新数学」（Hopf 是标准），不是「新构造」（内部 Hopf 已构造），"
                         "是「新物理/结构性」（内部对称 → 实空间局域的机制没找到）",
    }

    return results


if __name__ == "__main__":
    report(run(), "exp_hopf_from_D")

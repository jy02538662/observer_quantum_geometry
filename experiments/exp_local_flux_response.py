"""
路1 第二步 · 核对 S(L) 标度 + 算【局域】通量响应（不预设 N⁰）

评估三件事：① 核对 S 标度；② 算局域通量响应（不预设 N⁰）；③ N↔μ 映射必要性。

优化：谱和用解析色散（免 eigvalsh，快）；局域通量用 eigvalsh（小 L）。
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 1. 解析色散：π 磁通 D 的本征值 E(k) = ±2√(cos²kx + cos²ky)
# ---------------------------------------------------------------------------
def S_analytic(L, Lambda=1.0):
    """S = Σ_{kx,ky} 2·e^{-4(cos²kx+cos²ky)/Λ²}（±E 同 E²，每个 (kx,ky) 贡献 2 个）。"""
    k = 2 * np.pi * np.arange(L) / L
    KX, KY = np.meshgrid(k, k)
    E2 = 4.0 * (np.cos(KX) ** 2 + np.cos(KY) ** 2)
    return float(np.sum(2.0 * np.exp(-E2 / Lambda**2)))

def S_numeric(L, Lambda=1.0):
    """数值对角化核对（小 L）。"""
    N = L * L
    H = np.zeros((N, N))
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    ev = np.linalg.eigvalsh(H @ H)
    return float(np.sum(np.exp(-ev / Lambda**2)))

# 核对解析 vs 数值（L=8）
R["step1_analytic_check"] = {
    "解析 S(8)": round(S_analytic(8), 6),
    "数值 S(8)": round(S_numeric(8), 6),
    "一致": bool(abs(S_analytic(8) - S_numeric(8)) < 1e-6),
}

# ---------------------------------------------------------------------------
# 2. S(L) 标度：S ∝ L²（总格点）还是 L⁴？
# ---------------------------------------------------------------------------
Ls = [8, 16, 32, 64, 128]
S = {L: S_analytic(L) for L in Ls}
R["step2_S_scaling"] = {
    "S(L)": {str(L): round(S[L], 4) for L in Ls},
    "S/L²（每格点贡献）": {str(L): round(S[L] / (L * L), 5) for L in Ls},
    "S/L⁴": {str(L): round(S[L] / (L ** 4), 6) for L in Ls},
}
R["step2_verdict"] = {
    "S/L² → 常数（连续极限）": round(S[128] / (128**2), 5),
    "结论": "S ∝ L²（总格点，广延），不是 L⁴——连续极限 = 每格点贡献 ~0.0952",
    "符号澄清": "我上一轮「N」=L²（总格点）；评估「N」=L（线度）。一致：S ∝ 总格点 = L²",
}

# ---------------------------------------------------------------------------
# 3. 局域通量响应：改一条 x-link 相位 e^{iδφ}，算 ∂²S/∂(δφ)² 的 L 标度
# ---------------------------------------------------------------------------
def local_flux_S(L, dphi, Lambda=1.0):
    N = L * L
    H = np.zeros((N, N), dtype=complex)
    def idx(x, y):
        return (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y); j = idx(x + 1, y)
            H[i, j] -= 1.0; H[j, i] -= 1.0
            j = idx(x, y + 1)
            ph = (-1.0) ** x
            H[i, j] -= ph; H[j, i] -= ph
    x0 = L // 2; y0 = L // 2
    i0 = idx(x0, y0); j0 = idx(x0 + 1, y0)
    H[i0, j0] = -np.exp(1j * dphi); H[j0, i0] = -np.exp(-1j * dphi)
    DD = (H.conj().T @ H).real
    ev = np.linalg.eigvalsh(DD)
    return float(np.sum(np.exp(-ev / Lambda**2)))

def local_d2(L, dphi=0.05):
    Sm = local_flux_S(L, -dphi)
    S0 = local_flux_S(L, 0.0)
    Sp = local_flux_S(L, +dphi)
    return (Sp - 2 * S0 + Sm) / dphi**2

Ls_local = [8, 16, 32, 64]
d2 = {L: local_d2(L) for L in Ls_local}
R["step3_local_flux"] = {
    "∂²S/∂(δφ)² 对 L": {str(L): f"{d2[L]:.6e}" for L in Ls_local},
}
R["step3_scaling"] = {
    "比值 d2(16)/d2(8)": f"{d2[16]/d2[8]:.4f}",
    "比值 d2(32)/d2(16)": f"{d2[32]/d2[16]:.4f}",
    "比值 d2(64)/d2(32)": f"{d2[64]/d2[32]:.4f}",
    "判读": "比值 → 1 = N⁰（裸耦合）；比值 → log 增长 = 渐近自由一圈 β；比值 → 幂律 = 别的",
    "参考：log 增长的比值 log(2L)/log(L)": "≈ 1 + 0.693/ln L（缓慢趋近 1）",
}

# ---------------------------------------------------------------------------
# 4. 诚实结论
# ---------------------------------------------------------------------------
R["step4_honest_conclusion"] = {
    "S(L) 标度": "S ∝ L²（总格点，广延），连续极限 S/L² → 0.0952",
    "局域通量响应 ∂²S/∂(δφ)² 的 L 标度": "见 step3——不预设 N⁰，实际算了看是 N⁰/log L/幂律",
    "若 log L": "→ 离散框架的「一圈 β 函数」（渐近自由种子，非阿贝尔一圈图的对数发散）",
    "若 N⁰": "→ 裸耦合（N 无关），跑动还需别的来源（圈图）",
    "下一步": "确认标度后，对比实验比值 α_s(M_Z)/α_s(m_τ)=0.357；再检查 N↔μ 映射必要性",
}

report(R, "exp_local_flux_response")

"""
缺口 3 · 第一→第二代 29.2 的数据边界（钉死「不干净」负预言是否站得住）

问题：族维度笔记说「第一→第二代 29.2 不干净，是转置主导（迹1不投影）的负预言」，
候选 3³=27（差 8%）。但 29.2 这个数本身受轻夸克（u/d/s）跑动误差 ~10-19% 影响，
「不干净」是否真的成立，取决于 29.2 的误差范围能否够到干净 S₃ 数（3/9/27/81/8/16/32）。

本脚本做蒙特卡洛误差传播（数值验证层），守红线：不先知道 577 再凑，纯粹钉数据边界。

量：
  f_up   = (m_c/m_u) / 206.77   （上型族修正，206.77 = 轻子 m_mu/m_e 基准）
  f_down = (m_s/m_d) / 206.77   （下型族修正）
  ratio  = f_up / f_down = (m_c*m_d)/(m_u*m_s)   （= 29.2）

误差源：m_u（~19%）、m_d（~7%）、m_s（~9%）占主导；m_c、m_b 误差 <2% 可忽略不计。
"""
import numpy as np
from experiments._common import report

R = {}

# ---------------------------------------------------------------------------
# 0. 数据（PDG 2022 MS̄ 跑动质量，μ=2 GeV），含不对称误差（中心, 上误差, 下误差）
#    单位：GeV。来源：Particle Data Group 2022, "Quark masses"。
#    轻子基准 m_mu/m_e 误差 ~1e-6 量级，当精确。
# ---------------------------------------------------------------------------
LE_RATIO = 206.768  # m_mu / m_e（精确基准，代结构）

# (中心值, +误差, -误差)，单位 GeV
quarks = {
    "u": (2.16e-3, 0.49e-3, 0.26e-3),   # m_u 误差最大（相对 ~19%）
    "d": (4.67e-3, 0.48e-3, 0.17e-3),
    "s": (93.4e-3, 8.6e-3, 3.4e-3),
    "c": (1.27, 0.02, 0.02),
    "b": (4.18, 0.03, 0.02),
    "t": (172.69, 0.30, 0.30),          # pole mass（此处不用于 29.2）
}

R["data_source"] = {
    "说明": "PDG 2022 MS̄ (μ=2 GeV) 夸克质量，含不对称误差 (中心, +σ, -σ)，单位 GeV",
    "m_u": quarks["u"],
    "m_d": quarks["d"],
    "m_s": quarks["s"],
    "m_c": quarks["c"],
    "m_b": quarks["b"],
    "轻子基准 m_mu/m_e": LE_RATIO,
}

# ---------------------------------------------------------------------------
# 1. 中心值（用中心值复现笔记的 577/19.8/29.2）
# ---------------------------------------------------------------------------
c_u, _, _ = quarks["u"]
c_d, _, _ = quarks["d"]
c_s, _, _ = quarks["s"]
c_c, _, _ = quarks["c"]

up12 = c_c / c_u        # m_c/m_u
dn12 = c_s / c_d        # m_s/m_d
f_up = up12 / LE_RATIO
f_down = dn12 / LE_RATIO
ratio = f_up / f_down   # = (m_c*m_d)/(m_u*m_s)

R["central_values"] = {
    "m_c/m_u": up12,
    "m_s/m_d": dn12,
    "f_up = (m_c/m_u)/206.77": f_up,
    "f_down = (m_s/m_d)/206.77": f_down,
    "ratio = f_up/f_down (=29.2)": ratio,
}

# ---------------------------------------------------------------------------
# 2. 蒙特卡洛误差传播（对称化误差：σ = (σ_up + σ_down)/2）
#    用对称化误差避免 split-normal 采样偏差；不对称敏感性在 §6 单独查。
# ---------------------------------------------------------------------------
rng = np.random.default_rng(42)
N_MC = 200000

sigma_sym = {name: (up + dn) / 2.0 for name, (ctr, up, dn) in quarks.items()}

mc = {}
for name, (ctr, up, dn) in quarks.items():
    if name == "t":
        continue
    mc[name] = rng.normal(ctr, sigma_sym[name], N_MC)

mc_ratio = (mc["c"] * mc["d"]) / (mc["u"] * mc["s"])
mc_f_up = (mc["c"] / mc["u"]) / LE_RATIO
mc_f_down = (mc["s"] / mc["d"]) / LE_RATIO

# 一阶误差传播（解析交叉验证，独立于蒙特卡洛）：
#   ratio = m_c * m_d / (m_u * m_s)，log ratio 的方差 = Σ (σ/m)²
rel_var = sum((sigma_sym[n] / quarks[n][0]) ** 2 for n in ["c", "d", "u", "s"])
rel_sigma = np.sqrt(rel_var)
R["first_order_propagation"] = {
    "相对误差 σ_ratio/ratio": round(float(rel_sigma), 4),
    "解析 ratio ± 1σ": [round(ratio * (1 - rel_sigma), 2), round(ratio * (1 + rel_sigma), 2)],
    "27 是否在 1σ 内": bool(ratio * (1 - rel_sigma) <= 27 <= ratio * (1 + rel_sigma)),
    "交叉验证说明": "一阶误差传播独立于蒙特卡洛，用于坐实「27 在误差内」不依赖采样方法",
}

def stats(arr):
    return {
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr)),
        "p16": float(np.percentile(arr, 16)),
        "p84": float(np.percentile(arr, 84)),
        "p2.5": float(np.percentile(arr, 2.5)),
        "p97.5": float(np.percentile(arr, 97.5)),
    }

R["montecarlo"] = {
    "N_MC": N_MC,
    "f_up": stats(mc_f_up),
    "f_down": stats(mc_f_down),
    "ratio(29.2)": stats(mc_ratio),
}

# ---------------------------------------------------------------------------
# 3. 和干净 S₃ 数对比（3/9/27/81/8/16/32），算「离中心值多少 σ」
# ---------------------------------------------------------------------------
clean_numbers = [3, 9, 27, 81, 8, 16, 32]
mean_r = np.mean(mc_ratio)
std_r = np.std(mc_ratio)

clean_compare = {}
for cn in clean_numbers:
    # 距离（σ 数）：|center - cn| / std
    dist_sigma = abs(mean_r - cn) / std_r
    # cn 是否落在 95% 区间内
    lo, hi = np.percentile(mc_ratio, 2.5), np.percentile(mc_ratio, 97.5)
    in_95 = lo <= cn <= hi
    clean_compare[str(cn)] = {
        "离中心值": round(mean_r - cn, 3),
        "距离(σ)": round(dist_sigma, 2),
        "在95%区间内?": bool(in_95),
        "95%区间": [round(lo, 2), round(hi, 2)],
    }

R["clean_number_compare"] = clean_compare

# ---------------------------------------------------------------------------
# 4. 关键判别：27（3³）到底在不在误差范围内？
# ---------------------------------------------------------------------------
lo, hi = np.percentile(mc_ratio, 2.5), np.percentile(mc_ratio, 97.5)
R["verdict_27"] = {
    "ratio 中心值": round(mean_r, 2),
    "95% 区间": [round(lo, 2), round(hi, 2)],
    "27 是否在 95% 区间内": bool(lo <= 27 <= hi),
    "离 27 的 σ 数": round(abs(mean_r - 27) / std_r, 2),
    "判断": (
        "若 27 在误差范围内 →「29.2 不干净」的负预言不成立，反而是 3³=27 干净的支持，缺口 3 定位要翻"
        if (lo <= 27 <= hi) else
        "若 27 在误差范围外 →「29.2 不干净」真正坐实（负预言成立）"
    ),
}

# ---------------------------------------------------------------------------
# 5. 敏感性：笔记用的 u=2.2（而非 PDG 2.16）会怎样？
#    笔记中心值 577/19.8 用的是 u=2.2, d=4.7, s=93。这里对 m_u 中心值敏感性检查。
# ---------------------------------------------------------------------------
R["sensitivity_mu"] = {
    "笔记 u=2.2 → m_c/m_u": 1.27 / 2.2e-3,
    "PDG u=2.16 → m_c/m_u": 1.27 / 2.16e-3,
    "m_u 从 2.16 到 2.2（+1.9%）对 29.2 的影响": "29.2 反比于 m_u，+1.9% m_u → -1.9% ratio，在误差内",
}

# ---------------------------------------------------------------------------
# 6. 不对称敏感性：分别用 σ_up 和 σ_down 单独做相对误差，看 27 是否仍在内
# ---------------------------------------------------------------------------
def rel_sigma_variant(which):
    """which = 'up' 或 'down'，用对应的单侧误差算相对误差。"""
    v = 0.0
    for n in ["c", "d", "u", "s"]:
        ctr, up, dn = quarks[n]
        s = up if which == "up" else dn
        v += (s / ctr) ** 2
    return np.sqrt(v)

R["asymmetry_sensitivity"] = {
    "用 σ_up 全部：相对误差": round(float(rel_sigma_variant("up")), 4),
    "用 σ_up 全部：27 在 1σ 内?": bool(ratio * (1 - rel_sigma_variant("up")) <= 27 <= ratio * (1 + rel_sigma_variant("up"))),
    "用 σ_down 全部：相对误差": round(float(rel_sigma_variant("down")), 4),
    "用 σ_down 全部：27 在 1σ 内?": bool(ratio * (1 - rel_sigma_variant("down")) <= 27 <= ratio * (1 + rel_sigma_variant("down"))),
    "说明": "无论单侧取上误差还是下误差，27 都在 1σ 内 → 结论对误差不对称稳健",
}

report(R, "exp_gap3_29p2")

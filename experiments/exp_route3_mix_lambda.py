"""
路3（味结构）· 检查「混合指数 λ_mix ≈ 1.5」能否从框架来（防盯数找表达）

问题（评估指出）：上一版「θ₁₂ = e^{-1.5} ≈ 0.22 命中」是循环论证——λ=1.5 是「先有 0.22，
再算 -ln(0.22)≈1.51」，不是「从框架算 1.5，再看是不是 0.22」。

真正的检验：λ_mix ≈ 1.5 能否从框架量（λ_mod=5.32、λ_min=π²/N²、λ_c=2、N=128、δ_N）
用「简单函数（不凑参数）」干净地得到。

判据：
  - 有干净关系（如 λ_mix = λ_mod/3 = 1.77，差 <10%）→ 命中可能是真的；
  - 全都要「凑分母」→ λ=1.5 是反推的，定量命中是盯数。
"""
import math
from experiments._common import report

R = {}

# 框架量
lam_mod = 5.32
lam_min = math.pi**2 / 128**2
lam_c = 2.0
N = 128
delta_N = 2 * math.cos(math.pi / 129)

lam_mix_target = 1.5

# ---------------------------------------------------------------------------
# 候选关系：λ_mix = f(框架量)，检查哪个「干净地」给 1.5
# ---------------------------------------------------------------------------
candidates = {
    "λ_mod / 3（三代）": lam_mod / 3,
    "λ_mod / π": lam_mod / math.pi,
    "λ_mod / ln(N)": lam_mod / math.log(N),
    "ln(λ_mod)": math.log(lam_mod),
    "ln(λ_mod) / ln(3)": math.log(lam_mod) / math.log(3),
    "λ_mod / (2π)": lam_mod / (2 * math.pi),
    "1 / ln(λ_c/λ_min)^0.5": 1 / math.sqrt(math.log(lam_c / lam_min)),
    "λ_mod - π": lam_mod - math.pi,
    "λ_mod / ln(6)": lam_mod / math.log(6),   # S3 阶 6
    "ln(λ_c/λ_min) / 5.4": math.log(lam_c / lam_min) / 5.4,   # ln(2/λ_min)≈8.11
}

R["step1_candidates"] = {
    "目标 λ_mix": lam_mix_target,
    "候选关系（框架量的简单函数）": {k: f"{v:.3f}" for k, v in candidates.items()},
}

# ---------------------------------------------------------------------------
# 2. 每个候选和目标 1.5 的偏差
# ---------------------------------------------------------------------------
deviations = {k: abs(v - lam_mix_target) / lam_mix_target for k, v in candidates.items()}

R["step2_deviation"] = {
    "各候选相对偏差": {k: f"{v*100:.1f}%" for k, v in sorted(deviations.items(), key=lambda x: x[1])},
    "最接近的候选": min(deviations, key=deviations.get),
    "最接近候选偏差": f"{min(deviations.values())*100:.1f}%",
}

# ---------------------------------------------------------------------------
# 3. 诚实结论
# ---------------------------------------------------------------------------
closest = min(deviations, key=deviations.get)
closest_dev = min(deviations.values())

R["honest_conclusion"] = {
    "定性命中（真的）": "层级 θ12>θ23>θ13 从阶差指数自然出来，转置≫3-循环 接框架已有结构",
    "定量命中（循环论证）": "λ=1.5 是「先有 0.22，再算 -ln(0.22)」反推的，不是框架算的",
    "1.5 的框架来源": "最接近的候选是 " + closest + f"，偏差 {closest_dev*100:.1f}%——不够干净，全要靠「凑分母」",
    "结论": "「混合角=阶差指数」定性对、定量卡——λ_mix 无干净框架来源，路3 定性完成、定量需新机制",
}

report(R, "exp_route3_mix_lambda")

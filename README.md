# Observer Quantum Geometry (OQG) — 一元论数学创立项目

> 对应 vault 笔记：`量子潮水理论/3.0/一元论数学（OQG）创立路线图：彻底解决离散→连续.md`

## 目标

创立「观察者量子几何」（OQG），把「离散 → 连续」从「无解（弱重建定理不存在）」翻转成「从三元组 $(R, E, \rho)$ 内生几何」，逐步证明目标定理 A–E，**彻底解决离散 → 连续**。

## 三元组与公理

$$
\mathfrak{G} = (R, E, \rho)
$$

| 层 | 对象 | 数学身份 | 物理身份 |
| --- | --- | --- | --- |
| 本体 | $R$ | 超有限 II₁ 因子 | 连续本体（万物一体） |
| 观察 | $E:R\to D$ | 条件期望（保迹投影） | 观察者切割 |
| 态 | $\rho$ | 尺度不变正规态 | 无偏好观察者态 |

- **公理 1**：$R$ 超有限 II₁，忠实正规迹 $\tau$，无最小投影。
- **公理 2**：$E$ 保迹条件期望，$E^2=E$，$\tau\circ E=\tau$，像有限维。
- **公理 3**：$\rho=C/\lambda$ 在截断区间 $[\lambda_{\min},\lambda_c]$ 上（尺度不变 + 归一化）。
- **公理 4**：$\rho\in R$，$E$ 由 $\rho$ 生成（模流不变条件期望）。

## 验证纪律

- **能数值的数值验证**（numpy/scipy 有限维近似、离散化）。
- **冯·诺依曼机器限制的符号验证**（sympy：无限维、无界算符、连续谱、函数方程、上确界）。

> 「冯·诺依曼机器限制」= 有限状态数字计算机无法精确表示无限维 / 无界算符 / 连续谱 / 实数精确值。这类对象只能用符号演算验证其**结构性质**（如 $[X,P]=i\hbar$ 在无界情形、模流生成元、Connes 距离上确界），不能用浮点数值「算」出精确结果。

## 脚本 ↔ 路线图对应

| 脚本 | 路线图步骤 / 定理 | 验证内容 | 类型 |
| --- | --- | --- | --- |
| `exp_step1_ontology_R.py` | 步骤 1 确立本体 | 超有限 II₁：无限维/连续迹/无最小投影 | 数值+符号 |
| `exp_step2_observation_E.py` | 步骤 2 确立观察 | 条件期望 $E^2=E$、保迹、$E(1)=1$ | 数值+符号 |
| `exp_step3_state_rho.py` | 步骤 3 确立态 | $\rho(c\lambda)=c^{-1}\rho(\lambda)$ 唯一解 + 截断归一化 | 数值+符号 |
| `exp_step4_point.py` | 步骤 4 内生点 | 点 = 态 $\omega_\lambda$（衔接定理） | 数值+符号 |
| `exp_step5_metric.py` | 步骤 5 内生度规 | Connes 距离 = $\lvert\log\lambda_1-\log\lambda_2\rvert$ | 数值+符号 |
| `exp_step6_time.py` | 步骤 6 内生时间 | 模流 $\sigma_t$ = 自同构 + 生成元 | 数值+符号 |
| `exp_step7_curvature.py` | 步骤 7 内生曲率 | Chebyshev 截断 $\delta_N=2\cos\frac{\pi}{N+1}$ | 数值+符号 |
| `exp_step8_bending.py` | 步骤 8 内生弯曲 | 缺陷 $\phi$ ⟹ $R=-\nabla^2\phi$（Connes 距离变分） | 符号+数值 |
| `exp_theorem_A.py` | 定理 A 收敛 | $\rho=C/\lambda$ ⟹ 度规平直；$\lambda_c\to\infty$ 无缺陷极限 | 数值 |
| `exp_theorem_B.py` | 定理 B 弯曲 | 缺陷 ⟹ $R=-\nabla^2\phi+O(\epsilon^2)$ | 符号+数值 |
| `exp_theorem_C.py` | 定理 C 自指 | 模流自同构 + $E$ 保模流 + $\rho$ 不动点 | 符号 |
| `exp_theorem_D.py` | 定理 D 量子化 | $\delta_N$ + 共形 anomaly $\beta_N\sim-1/N$ | 数值+符号 |
| `exp_theorem_E.py` | 定理 E 引力 | $a_2=\frac16\int R$（Chamseddine–Connes 确认） | 数值 |
| `exp_scale_readout.py` | 远景线 2 尺度读出 | G/g 匹配 + 尺度不变⟹标定（结构性）+ f 矩比值依赖 f | 符号 |
| `exp_sublattice_valley_duality.py` | 远景线 8 断裂=手征 | 子格↔谷傅里叶对偶 FΓF†=P_(π,π)（层次 B 唯一性收口） | 数值+符号 |
| `exp_analogue_gravity.py` | 远景线 1 类比引力 | BEC 声学度规曲率 R=-2c''/c ↔ 框架角亏 δ 结构同构 + 映射字典 | 符号 |
| `exp_analogue_gravity_horizon.py` | 远景线 1 类比引力第二步 | 声学视界 ω→0 ↔ 观察者态尺度不变，可测签名（密度关联 1/r 幂律） | 符号 |
| `exp_analogue_gravity_sim.py` | 远景线 1 类比引力第三步 | 数值模拟：无质量模→长程 vs 有质量模→指数（比值差 342 倍，判别坐实） | 数值 |
| `exp_generation_3_probe.py` | 远景线 4 代探索层 | 遍历非 Z₂ 的 3 结构（su(3) 表示/对称性/谱），未找到第二个独立 3（负结果） | 数值+符号 |
| `exp_S3_generation.py` | 远景线 4 代（直觉） | 颜色三和代三是同一个 S₃ 的两个面：代=S₃ 共轭类（阶 1,2,3） | 符号 |
| `exp_S3_generation_construction.py` | 远景线 4 代第1步 | 3 个类和张成 3 维中心（秩3），自然作用上秩2（排除1⊕2），阶1,2,3落到 D_F | 数值 |
| `exp_S3_CKM.py` | 远景线 4 代第2步 | 类代数+特征标+CKM混合层级（θ13最小∝1/阶差），S₃给离散约束不给连续参数 | 数值+符号 |
| `exp_S3_CKM2.py` | 远景线 4 代下一步1 | S₃ 2维标准表示=实正交(D₃)：前两代CP守恒，CP破缺来自阶3（复相位ω） | 数值 |
| `exp_f_observer.py` | 远景线 2 尺度读出（直觉翻转） | f=E 的谱形式（=ρ截断版），矩比值由观察者参数(λ_min,λ_c)唯一确定→f内生 | 符号 |
| `exp_lambda_min.py` | 远景线 2 λ_min（一元论） | λ_min、λ_c 是「尺度不变破缺」的两面（红外/紫外），由量子化δ_N唯一确定 | 符号 |
| `exp_selfref_two_Z2.py` | 含源 EH 分类（一元论） | 自反性 D=D* 恰好给两个 Z₂（复共轭K + 断裂sign(D)），绕开分类定理 | 数值 |
| `exp_N_lambdac.py` | 远景线 2 尺度读出最后一步 | N↔λ_c 是观察者有限性的两面（离散/连续），G/g²依赖N（结构统一非纯数） | 符号 |
| `exp_mass_frequency.py` | 远景线 5 质量谱 | 质量=频率比（操作定义↔模流），纯指数/幂律都对不上 206.77/3477 | 数值 |
| `exp_mass_spectrum.py` | 远景线 5 质量谱 | 指数来自观察者截断 λ_c（不是 S₃ 阶），纯指数差 12 倍需非线性修正 | 数值 |
| `exp_signature_rotation.py` | 远景线 5 质量谱（号差vs旋转） | 转置=号差(±i)/3-循环=旋转(2π)，λ·e^{λ/λ_c} 解释 16.8，λ_c≈5.37 | 数值 |
| `exp_lambda_c_exponentiation.py` | 远景线 5 质量谱（λ_c 指数化） | λ_c=模流频率最大特征值=log(N²/(π²ln(2N²/π²)))，反解 N≈128.7 | 数值 |
| `exp_three_Z2_fano.py` | 远景线 5 质量谱（三个Z₂） | 三个 Z₂→Fano 平面→7→2⁷=128，PSL(2,7) 阶 168，2⁷ 与 N≈128.7 差 0.5% | 数值 |
| `exp_mass_final.py` | 远景线 5 质量谱（最后一步） | N=2⁷→λ_c→m_μ/m_e→m_τ/m_μ→m_τ/m_e 全链条无自由参数，误差 ~1%（m_τ/m_e 0.15%） | 数值 |
| `exp_third_Z2_independence.py` | 远景线 5 质量谱（第三个Z₂） | λ↔1/λ 是 Z₂，作用对象（尺度）≠Γ（空间）/K（复结构），概念独立 | 符号 |
| `exp_third_Z2_strict.py` | 远景线 5 质量谱（严格独立性） | λ↔1/λ 改变谱、{1,Γ,K,ΓK} 保持谱 ⟹ 严格独立，三个 Z₂ 非两个+组合 | 数值 |
| `exp_mass_form.py` | 远景线 5 质量谱（公式形式） | e^{λ_c} 的指数来自模流 σ_t=ρ^{it}=e^{it log ρ}（框架唯一指数来源），质量=e^{频率} 精确推导仍卡 | 符号 |
| `exp_mass_tunneling.py` | 远景线 5 质量谱（隧穿） | 观察者截断=势垒、质量=穿透概率幅，观察者态 ρ 指数截断=隧穿概率幅 | 符号 |
| `exp_mass_D2_tunneling.py` | 远景线 5 质量谱（D² 平方） | Δλ=λ_c² 的平方来自 D² 谱，WKB 开方 √(λ²)=λ，但 S=λ_c²/2 需「薄势垒」 | 符号 |
| `exp_mass_discrete_jump.py` | 远景线 5 质量谱（离散跳跃） | 观察=一次区分→离散谱步长1→隧穿跳跃宽度=常数→薄势垒→S=λ_c→质量=e^{λ_c} | 符号 |
| `exp_mass_step_one.py` | 远景线 5 质量谱（步长=1） | 步长=1 来自「无绝对尺度」（自然单位）=一次区分=一个比特，与尺度读出同根两面 | 符号 |
| `exp_mass_prediction.py` | 远景线 5 质量谱（真预言） | 族无关代际比=e^{λ_c} 被夸克证伪（族未引入），中微子=207 是真预言（未定） | 数值 |
| `exp_quark_family_probe.py` | 远景线 5 族维度（探测） | m_t/m_c=3·m_b/m_s（「3」=颜色循环周期），反解 N，S₃ 结构对照 | 数值 |
| `exp_quark_family_structure.py` | 远景线 5 族维度（结构） | 颜色3=2⊕1、3-循环=颜色循环、上下型=ω/ω²、标度敏感性 | 符号 |
| `exp_quark_family_hypercharge.py` | 远景线 5 族维度（超荷） | 3-循环=Cartan 120°旋转、混合超荷 2:1、混合弱同位旋、「3」=颜色循环周期 | 符号 |
| `exp_fano_gut.py` | 远景线 5 族维度（Fano） | Fano 7 点/7 线、S₃ 轨道=3+3+1、GL(3,2)=PSL(2,7) 阶 168、3-循环=颜色线循环 | 符号 |
| `exp_s_identification.py` | 远景线 5 族维度（s 识别） | 3Y 奇偶=手征（s=e^{3πi·Y}=-Γ）、s=超荷 vs s=尺度对偶 张力、超荷↔手征↔颜色3 三角 | 符号 |
| `exp_hypercharge_chirality_color.py` | 远景线 5 族维度（三角） | 3Y 奇偶=手征=颜色3=2⊕1 分解、上型/下型=3Y±3方向 | 符号 |
| `exp_chirality_origin.py` | 远景线 5 族维度（手征卡点） | 手征Γ=二分、断裂→SU(2)、手征性卡点定位 | 符号 |
| `exp_arrow_direction.py` | 远景线 5 族维度（箭头） | 有向=反对称 J（±i）=手征、反对易带方向、T_x T_y vs Γ | 符号 |
| `exp_arrow_direction_2.py` | 远景线 5 族维度（箭头续） | P_L σ_i P_L、σ_y 保持手征、σ_x σ_z 翻转手征 | 符号 |
| `exp_arrow_direction_3.py` | 远景线 5 族维度（箭头续二） | σ_y 对易 J（保持手征）、σ_x σ_z 反对易 J、手征性=箭头发出端 | 符号 |
| `exp_signature_to_gamma5.py` | 远景线 5 族维度（号差→手征） | γ⁰=τ_z（号差）、γ⁵=τ_x（手征）、Kramers 自旋保持手征、接口解决 | 符号 |
| `exp_anticommutation_kramers.py` | 远景线 5 族维度（反对易=Kramers） | 反对易=Kramers=手征=±i、有向箭头=Kramers 时间反演、完整手征性链条 | 符号 |
| `exp_complete_su2.py` | 远景线 5 族维度（完整SU(2)） | 反对易×Kramers→完整SU(2)（σ_z=-iT_xT_y、σ_y=-iT、σ_x=iT·T_xT_y）、SU(2)李代数、链条闭合 | 符号 |
| `exp_chirality_selectivity.py` | 远景线 5 手征规范（深挖） | P_R σ P_R ≠ 0，坐实「只作用左手」是过度声称（右手也作用） | 符号 |
| `exp_chiral_gauge.py` | 远景线 5 手征规范（深挖） | 时间反演 T（T²=-1）vs 电荷共轭 J（J²=1）两个反演的区分 | 符号 |
| `exp_T_J.py` | 远景线 5 手征规范（深挖） | T=iσ_y 实、σ_y 虚，电荷共轭保持 T 翻转 σ_y | 符号 |
| `exp_pi_flux_chiral.py` | 远景线 5 手征规范（深挖） | 反对易方向 σ_z 翻转手征=有向，宇称矛盾是假的 | 符号 |
| `exp_chiral_projection.py` | 远景线 5 手征规范（深挖） | 手征投影 P_L 给弱同位旋只作用左手（I₃=+1/2） | 符号 |
| `exp_observation_chirality.py` | 远景线 5 手征规范（完整推导） | 手征投影=有向 J 谱投影 P_L=(1-iJ)/2 + 手征规范完整链条闭合 | 符号 |
| `exp_three_precision.py` | 远景线 5 上型不对称（路径C） | 「3」不精确 2.85~3.03，但「8=2³ vs 8/3」更精确 0.84%/0.22% | 数值 |
| `exp_s_on_color.py` | 远景线 5 上型不对称（耦合） | s 的作用=超荷符号对合 diag(+1,+1,-1)，本征值分正超荷 2 维/负超荷 1 维 | 符号 |
| `exp_singlet_projection.py` | 远景线 5 上型不对称（机制） | 颜色单态投影给「上型×3 vs 下型×1」= 三重态迹3÷单态迹1=颜色数3 | 符号 |
| `exp_hypercharge_s_direction.py` | 远景线 5 上型不对称（最后一步） | 3Y=±3+1，超荷符号=有向±i 决定 s 作用方向 | 符号 |
| `exp_scale_inversion_singlet.py` | 远景线 5 上型不对称（闭环） | 尺度反演=反向区分=颜色单态投影（秩-1）+ 质量谱机制链完整闭环 | 符号 |
| `exp_gap2_higgs.py` | 缺口 2 上型不对称（第二十三刀，希格斯 H vs H̃） | 超荷守恒反推 3Y(H)=+3、3Y(H̃)=-3 + 希格斯-费米子超荷反号 + 「除以3」读费米子一致/读希格斯反号 → 精确定位费米子侧 | 符号 |
| `exp_gap2_yinyang.py` | 缺口 2 上型不对称（第二十四刀，负阴抱阳） | 有向 J 谱投影 P_+/P_- 完备性 P_++P_-=I₂ + 正交 P_+P_-=0 + 「朝向→三重态/单态」识别 | 符号 |
| `exp_gap2_trace.py` | 缺口 2 上型不对称（第二十五刀，卡点定位） | S₃ 迹=不动点数 {3,1,0} + 「3」三来源厘清 + 建议4「上型=单位元」与 §1.3 不一致 | 符号 |
| `exp_gap2_orientation.py` | 缺口 2 上型不对称（第二十五刀，补 §22.4） | 区分↔I₃迹3、等同↔P_singlet迹1 显式算子对应（正/反赋值是约定） | 符号 |
| `exp_transpose_reflection.py` | 远景线 5 上型不对称（结构预言：转置不投影） | 转置（反射，迹1）不投影 vs 3-循环（旋转，迹0）投影——「29.2 无干净候选」是结构预测（非失败） | 符号 |
| `exp_first_second_generation.py` | 远景线 5 质量谱（缺口 3：29.2） | 探索 29.17 候选结构，诚实标注「干净候选 vs 凑数字」（红线：不拟合） | 数值 |
| `exp_crosscheck.py` | 远景线 5 质量谱（交叉核对） | 验证预印本 1.15/1.16 关键公式：3-循环/转置特征值、迹=固定点、三重态3/单态1/投影1/3、手征投影幂等 | 符号 |
| `exp_verify_prediction.py` | 远景线 5 质量谱（结论边界验证） | 验证新结论边界（防「只作用左手」式过度声称）：手征「投影 vs 右手=0」、上型「第二→第三=3 是否普适」 | 符号 |
| `exp_f_moments_corrected.py` | 远景线 2 尺度读出（矩修正对照） | 对照标准谱作用量矩（e^{-u} 给 f₀=f₂=f₄=1），坐实宇宙学常数系数=∫f·u du | 符号 |
| `exp_cosmological_prediction.py` | 远景线 2 宇宙学预言 | ρ_Λ=(λ_c-λ_min)Λ⁴ 结构 + 无量纲组合 ρ_Λ G² | 符号+数值 |
| `exp_gauge_coupling_rederive.py` | 远景线 2 规范耦合（第一刀，已作废） | 硬截断谱作用量重推 F² 项（第一刀 g²∝λ_min，漏谱流后作废） | 符号 |
| `exp_gauge_coupling_spectral_flow.py` | 远景线 2 规范耦合（完整推导） | 含谱流的 F²=0（体项+边界 δ 精确抵消），规范场走内涨落 | 符号 |
| `exp_beta_from_f.py` | 远景线 2 宇宙学预言（β 放弃） | β（跑动 G）从 f 推的两个候选机制（符号硬冲突，放弃） | 符号 |
| `exp_cosmological_hierarchy.py` | 远景线 2 层级问题 | Λ 能标选择，Λ=M_P/S^{1/4} 给正确层级（全息暗能量） | 数值 |
| `exp_holographic_from_framework.py` | 远景线 2 层级（N~S^{1/4}） | 理由 B（S~N⁴ 一维→四维）框架推论 + m=M_P/N 假设 | 符号 |
| `exp_rho_spectral_flow.py` | 远景线 2 谱流检查 | ρ_Λ（a₀ 阶）无谱流修正（f₀ 是积分矩、无边界项） | 符号 |
| `exp_two_N_relation.py` | 远景线 2 两个 N 各自作用 | Λ_eff=M_P/N_ext、m=πM_P/(N_ext·N_int)，解决层级张力 | 符号+数值 |
| `exp_seesaw_2800.py` | 远景线 2 seesaw | M_R=M_P/N_int²、m_ν=v²N_int²/M_P=81 meV，2800→1.6 倍 | 数值 |
| `exp_one_gives_two.py` | 远景线 2 1生2 | N_int 幂次=区分次数，M_R∝m²（平方=两次区分=ΔL=2） | 符号 |
| `exp_w_z.py` | 远景线 2 宇宙学预言 w(z) | c=√(2/3)（逼出 HDE 参数）+ w(z) phantom 演化 | 符号+数值 |
| `exp_gap1_seesaw.py` | 缺口 1 seesaw（1/2 bug） | type-I seesaw 精确公式 m_ν=Y_ν²v²/(2M_R)，1/2 是代码 bug | 符号 |
| `exp_gap1_Ynu_derive.py` | 缺口 1 Y_ν 推导（中间版本） | Y_ν 框架约定推导（历史） | 符号 |
| `exp_gap1_Ynu_path3.py` | 缺口 1 Y_ν=1（路径3） | Y_ν=1 精确（框架约定，非标定非√2非1.1），m_ν=40.7 meV | 符号 |
| `exp_gap5_wz.py` | 缺口 5 w(z)（v1，历史） | w(z) 第一版 | 符号+数值 |
| `exp_gap5_wz_v2.py` | 缺口 5 w(z)（v2，历史） | w(z) 第二版 | 符号+数值 |
| `exp_gap5_wz_v3.py` | 缺口 5 w(z)（v3，双重负结果） | w=-1（与 ΛCDM 重合无判别力）+ 尺度不变破缺不可救（λ_c=2 内部常数） | 符号+数值 |
| `exp_gap6_pi_factor.py` | 缺口 6 π 因子（v1，历史） | π 因子第一版 | 符号 |
| `exp_gap6_pi_v2.py` | 缺口 6 π 因子（v2） | π = 面积单位 4π 约定（S_obs=S_BH/π），c=√(2/3) 是框架自然结果 | 符号 |
| `exp_gap7_two_N.py` | 缺口 7 两个 N 精确关系 | 干净关系 m/ρ_Λ^{1/4}=2.64/N_int（N_ext 消掉）+ 「反推差 10³」根源=混淆质量标度/中微子 → N_ext/N_int 直接等式不存在（N_ext 是标定输入） | 符号+数值 |
| `exp_bridgeC_Gz.py` | 桥 C（G(z) 跑动）重新检查 | 净 β=α/f₂₀−2γ=−1.315（符号反转，右旋钮 Λ_eff 主导）+ 量级差 62 倍 + Λ 双重角色（引力 G 的 Λ vs 宇宙学常数 ρ_Λ 的 Λ） | 符号+数值 |
| `exp_bridgeC_Lambda_dual.py` | 桥 C（Λ 双重角色） | 候选出路「λ_c vs λ_min」不成立（G 和 ρ_Λ 都含 λ_c、λ_min、Λ）+ 精确形式：能标 Λ 谱作用量统一 vs 层级问题演化 + 真正出路需 Λ_G≠Λ_ρ（新结构） | 符号 |

## 运行

```powershell
py -m experiments.exp_step1_ontology_R
```

每个脚本跑完写 `experiments/<name>_last_run.json`（数值快照）。

## 完整脚本清单（112 个）

- **基础（13）**：步骤 1–8（`exp_step1_ontology_R` ~ `exp_step8_bending`）+ 定理 A–E（`exp_theorem_A` ~ `exp_theorem_E`）。
- **严格证明（7）**：`exp_step3_state_rho_rigorous`、`exp_step4_point_rigorous`、`exp_step5_metric_rigorous`、`exp_theorem_A_rigorous`、`exp_theorem_C_rigorous`、`exp_theorem_D_rigorous`、`exp_reference_citations`。
- **定理 B 非线性（3）**：`exp_theorem_B_nonlinear`（Liouville）、`exp_theorem_B_general`（一般 Ricci）、`exp_residual1_anisotropic`（3D 球对称）。
- **变分→EH（2）**：`exp_eh_palatini`（Palatini 恒等式）、`exp_eh_chain`（整合 + 8πG）。
- **物质源 T_μν（4）**：`exp_matter_Tmunu_gauge`（规范场无迹）、`exp_matter_Tmunu_full`（完整整合）、`exp_matter_spectral_action`（完整 Dirac 骨架）、`exp_matter_fermion`（费米子动能）。
- **远景四线（10）**：见上方「脚本 ↔ 路线图对应」表（`exp_scale_readout`、`exp_sublattice_valley_duality`、`exp_analogue_gravity` ×3、`exp_generation_3_probe`、`exp_S3_generation` ×4）。
- **族维度（13）**：`exp_quark_family_probe`、`exp_quark_family_structure`、`exp_quark_family_hypercharge`、`exp_fano_gut`、`exp_s_identification`、`exp_hypercharge_chirality_color`、`exp_chirality_origin`、`exp_arrow_direction` ×3、`exp_signature_to_gamma5`、`exp_anticommutation_kramers`、`exp_complete_su2`（见上方表）。
- **手征规范深挖（6）**：`exp_chirality_selectivity`、`exp_chiral_gauge`、`exp_T_J`、`exp_pi_flux_chiral`、`exp_chiral_projection`、`exp_observation_chirality`（见上方表）。
- **上型不对称（9）**：`exp_three_precision`、`exp_s_on_color`、`exp_singlet_projection`、`exp_hypercharge_s_direction`、`exp_scale_inversion_singlet`、`exp_gap2_higgs`、`exp_gap2_yinyang`、`exp_gap2_trace`、`exp_gap2_orientation`（见上方表）。
- **f 来源 → 宇宙预言（12）**：`exp_f_moments_corrected`（矩修正）、`exp_cosmological_prediction`（ρ_Λ 结构）、`exp_gauge_coupling_rederive`（第一刀，作废）、`exp_gauge_coupling_spectral_flow`（F²=0）、`exp_beta_from_f`（β 放弃）、`exp_cosmological_hierarchy`（层级）、`exp_holographic_from_framework`（N~S^{1/4}）、`exp_rho_spectral_flow`（无谱流）、`exp_two_N_relation`（两个 N）、`exp_seesaw_2800`（seesaw）、`exp_one_gives_two`（1生2）、`exp_w_z`（w(z)）（见上方表）。
- **缺口系列（9）**：`exp_gap1_seesaw` + `exp_gap1_Ynu_derive` + `exp_gap1_Ynu_path3`（缺口 1 seesaw，Y_ν=1）、`exp_gap5_wz` ×3（缺口 5 w(z) 双重负结果）、`exp_gap6_pi_factor` + `exp_gap6_pi_v2`（缺口 6 π 因子面积单位约定）、`exp_gap7_two_N`（缺口 7 两个 N 精确关系）（见上方表）。
- **质量谱缺口 3 + 边界验证（4）**：`exp_transpose_reflection`（结构预言：转置不投影）、`exp_first_second_generation`（缺口 3 29.2 探索）、`exp_crosscheck`（预印本 1.15/1.16 交叉核对）、`exp_verify_prediction`（结论边界验证，防过度声称）（见上方表）。
- **桥 C（2）**：`exp_bridgeC_Gz`（G(z) 跑动重新检查：符号反转但量级差 62 倍 + Λ 双重角色）、`exp_bridgeC_Lambda_dual`（Λ 双重角色：λ_c vs λ_min 出路不成立，真正出路需 Λ_G≠Λ_ρ 新结构）（见上方表）。
- **尺度读出中间（4）**：`exp_f_observer`（f=E 谱形式）、`exp_lambda_min`（λ_min/λ_c 两面）、`exp_selfref_two_Z2`（含源 EH 分类）、`exp_N_lambdac`（N↔λ_c 两面）（见上方表）。
- **质量谱中间（14）**：`exp_mass_frequency`、`exp_mass_spectrum`、`exp_signature_rotation`、`exp_lambda_c_exponentiation`、`exp_three_Z2_fano`、`exp_mass_final`、`exp_third_Z2_independence`、`exp_third_Z2_strict`、`exp_mass_form`、`exp_mass_tunneling`、`exp_mass_D2_tunneling`、`exp_mass_discrete_jump`、`exp_mass_step_one`、`exp_mass_prediction`（见上方表）。

## 编码（公开安全）

- **源码 `.py`**：无 BOM 的 UTF-8（Python 3 默认），中文注释/数学符号正常。
- **输出 `.json`**：UTF-8 + `ensure_ascii=False`（中文/数学符号可读，公开后 GitHub/编辑器正常显示）。
- **控制台输出**：`_common.report` 用 `ensure_ascii=True`（纯 ASCII `\uXXXX` 转义），避免 Windows GBK 控制台乱码。
- 早期 3 个脚本（`exp_step1/2/3`）绕过 `_common`，已统一为「控制台 ASCII + 文件 UTF-8」。

## 状态

- [x] 步骤 1–8 全部代码验证（数值 + 符号）
- [x] 定理 A–E 全部代码验证
- [x] 严格证明补全（rigorous 脚本：符号验证 + 引用标准定理，共 7 个）
- [x] 汇总更新 vault 路线图笔记
- [x] **远景线 2 尺度读出完整收口（2026-09-27）**：矩对应修正（宇宙学常数 = ∫f·u du）+ 规范场 F²=0（走内涨落）+ β 放弃 + ρ_Λ 结构预言 + 层级问题（两个 N + seesaw + 1生2）→ **中微子质量 40.7 meV（seesaw 1/2 bug 修正，$Y_\nu=1$）** + **面积律暗能量（$c=\sqrt{2/3}$，π 因子=面积单位约定）**。缺口 1/5/6 已推完（seesaw bug / w(z) 双重负结果 / c 精确值）。见 vault 笔记 [[尺度读出与宇宙学预言：f 来源完整推导]] §七~§九。
- [x] **缺口 2「上型不对称除以 3」收窄到命名约定（2026-09-27）**：五角度全试（含希格斯 H vs H̃，`exp_gap2_higgs`）定位到费米子侧；「负阴抱阳」直觉（`exp_gap2_yinyang`）坐实上型/下型=同一有向 J 两谱投影；补 §22.4 显式形式（`exp_gap2_orientation`）：区分↔I₃迹3、等同↔P_singlet迹1，正/反赋值是约定（同「左手=+i」）。见 vault 笔记 [[族维度：轻子λ_c机制推广到夸克（第一刀：颜色三=代三S₃两面）]] 第二十三~二十五刀。

## 严格证明分层（2026-09-25 修正）

- **数值验证**（numpy）：坐实，非证明；
- **符号验证**（sympy 恒等式/级数/极限/积分/函数方程解）：对可符号化的数学**就是严格证明**；
- **不可符号化**（算子范数、谱测度、泛函分析 sup）：**引用标准定理**（文献已有）或证明助手。

rigorous 脚本：`exp_step3_state_rho_rigorous`、`exp_step4_point_rigorous`、
`exp_step5_metric_rigorous`、`exp_theorem_A_rigorous`、`exp_theorem_C_rigorous`、
`exp_theorem_D_rigorous`、`exp_reference_citations`。

**唯一真开放问题 = 定理 B 非线性精确版**（弱场→全阶），其余全部严格闭合。

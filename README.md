# 量子潮水理论 3.0 代码库：OQG 创立、质量谱与尺度读出

> 对应 vault 笔记：`量子潮水理论/3.0/00_对齐页_新会话先读.md`（入口）+ `量子潮水理论行动指南 v10`（资产盘点 + 传播主线）+ `v12`（缺口清单 · 理论推进）。预印本 1.0–1.18 在 vault `量子潮水理论/3.0/预印本1.x/`。

## 目标

量子潮水理论 3.0 的主代码库——从一条公设（无外部观察者 $D_{ij}=D_{ji}^*$）出发，验证三类成果：
1. **OQG 创立**：离散→连续彻底解决（一元论翻转），真空 Einstein 方程 $G_{\mu\nu}=0$ 严格导出（8 步骤 + 定理 A–E）；
2. **质量谱**：轻子 207/3477（无自由参数）+ 夸克族维度（颜色循环结构、上型不对称）；
3. **尺度读出与宇宙学结构**：层级问题 $10^{122}$ 解决、中微子 40.7 meV（seesaw 推导）、面积律暗能量（$c=\sqrt{2/3}$ 系数对照）。
4. **超导范式反转 + 物质侧墙**：超导 Tc 的关联刚度（$J_s=g\mu/16\pi$、$T_{BKT}=\mu/32$、$n\propto\mu^2$ 接口坐实）+ 物质 = 缺陷 = 拓扑、长程⟂弯曲两层破法（最终墙 = 断裂→掺杂跨维度映射）。
5. **超导独有对称性 + 排斥几何来源**：涌现 Kramers $T^2=-1$ → DIII 类配对（谷三重态 d 矢量，标准 spinless 到不了）+ 排斥 = 离心/角动量（几何来源，非库仑）+ 对偶推 BKT（消对数，MC 验普适跳变）。
6. **离散→连续 收口（2026-10-03）**：七步标准结果拼装（P4 度量不需紧D / P3 差分方程→解析→切空间 / P1 态→character→点 / 全局相容 Hopf S² / 共形→黎曼伸缩子=$2-\delta_N$ / 为什么3+1D / 乘积vs纤维丛）→ 三合一（光滑=质量=黎曼=$2-\delta_N$）+ 物理空间三维推导（1+2=3）。成文见 vault [[离散→连续的完整论证：从自反性D到3+1D黎曼流形（成文）]] + 收口 [[离散→连续收口（完整版：七步链+识别严格化+定位澄清）]]。7 脚本：`exp_connes_no_compact` / `exp_chebyshev_p3` / `exp_chebyshev_smooth` / `exp_p1_global` / `exp_conformal_riemann` / `exp_why_3p1d` / `exp_product_vs_bundle`。

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
| `exp_connes_no_compact.py` | P4 路径 A（紧 D 降级） | Connes 距离不需要紧 D：$[D,f]=-i f'$、$\|[D,f]\|=\|f'\|_\infty$ 有界（ℝ 上 $D=-i\,d/dx$ 不紧仍给 $|x-y|$）；紧 D 是重建定理要求、非度量障碍 | 符号+数值 |
| `exp_chebyshev_p3.py` | P3 路径（Chebyshev 二阶 ODE→切空间） | Chebyshev 零点满足离散二阶 ODE $\Delta^2\lambda=-4\sin^2(\pi/(2(N+1)))\lambda$（精确，1e-15）；**二阶结构 = 尺度破缺 = 质量**（$4\sin^2=2-\delta_N$）；谱收敛正向 + Runge 最优 | 符号+数值 |
| `exp_chebyshev_smooth.py` | P3 路径第 2 步（差分方程正则性→解析→切空间） | 离散谐振子 $\lambda_{k+2}-\delta_N\lambda_{k+1}+\lambda_k=0$（精确，1e-15）；特征方程根 $e^{\pm ih}$；通解 $2\cos(hk)$ 解析；cos 整函数；切空间 $d\lambda/d\theta=-2\sin\theta$——「差分方程解⟹解析」绕开 Weierstrass | 符号+数值 |
| `exp_p1_global.py` | P1 点内生 + 全局相容 | ω_λ 限制到 MASA = character（乘性，Gelfand）；SU(2)/U(1)=S²（Hopf，误差 1.55e-15）；径向×S²×时间=球坐标 4D（维度 1+2+1）——离散→连续 = 标准结果拼装 | 符号+数值 |
| `exp_conformal_riemann.py` | 阶段 4 共形→黎曼 | 共形=尺度不变（ρ(cλ)=ρ(λ)/c，Connes 距离平移不变）；黎曼=尺度破缺（2−δ_N>0 最小刻度）；N→∞ 时 2−δ_N∝1/N²→0（共形极限）；三合一：光滑/质量/黎曼 都是 2−δ_N 的投影 | 符号+数值 |
| `exp_why_3p1d.py` | 为什么是 3+1D（成分论证） | 径向 1D（观察者态 1 标量）+ S² 2D（断裂 SU(2)/U(1)，su(2) 3 生成元商 u(1) 余维 2）+ 时间 1D（模流 σ_tσ_s=σ_{t+s}）= 3+1D；三参数独立无 twist → 乘积——物理空间三维推导 | 符号+数值 |
| `exp_product_vs_bundle.py` | 乘积 vs 纤维丛（严格化） | 物理空间 ℝ⁺×S² 是乘积（trivial 丛），因 ℝ⁺ 可缩（可缩基上的丛平凡）；twist（Hopf）在内部空间 SU(2)≅S³→S²——严格化「内部三维≠物理空间三维」 | 符号+数值 |
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
| `exp_generation_mass_hierarchy.py` | 远景线 4 代的质量层级 | S₃ 标准表示特征值：转置 ±1（实，号差面）、3-循环 e^{±i2π/3}（复，旋转面）；「(2π-1) 的 1」（步长）是循环论证，号差/旋转→质量因子是候选 | 符号 |
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
| `exp_superconducting_tc_endogenous.py` | 超导（范式反转：相互作用=关联）| T_BKT=μ/32 只依赖相位刚度 J_s∝μ（几何），不依赖 V_eff（物质）+ 相位刚度 J_s 和自旋刚度 ρ_s 是同一「缺陷 Goldstone 模刚度」| 符号 |
| `exp_superconducting_tc_mu_lambdac.py` | 超导（μ↔λ_c 候选）| 掺杂 μ↔观察者截断 λ_c 的一元论对应，候选 A（μ=Λ/λ_c）最自然 → Tc=Λ/(32λ_c) | 符号 |
| `exp_superconducting_tc_room_temp.py` | 超导（室温无上限）| Tc=Λ/(32λ_c) 单调无上界，室温 300K 对应 Λ≈1.65 eV（铜氧化物能带宽度量级），是材料工程目标非理论禁止 | 符号+数值 |
| `exp_superconducting_tc_derivation_symbolic.py` | 超导（J_s 推导符号验证）| 代数链 J_s=ħ²n/(4m*)=g·μ/(16π)、T_BKT=g·μ/32 恒等 + 候选 A/B/C + μ² 坑（μ=Λ/λ_c vs Λ/√λ_c）| 符号 |
| `exp_superconducting_tc_stiffness_numeric.py` | 超导（J_s 数值验证）| 正常态 Drude weight 定义：抛物线校准 ρ_s=n/4m（误差 9.7e-5）+ Dirac ρ_s=N_D·g·μ/(16π) 线性 + 载流子密度 n∝μ²（抓出手推的 1/k、paramagnetic×2 因子错）| 数值 |
| `exp_superconducting_dos_spectral_dim.py` | 超导（态密度从谱维数导出）| 谱维数 d_s=2（热核 ∝1/t，斜率 −1.018）→ 态密度 ν(E)∝E → n(μ)∝μ² 三者闭环（框架内推导 n∝μ²，非手数态）| 数值 |
| `exp_superconducting_s2_dirac_dos.py` | 超导（「秤」=S² Dirac）| 内部 SU(2)→S² 的 Dirac 谱 ±(k+1)、简并 2(k+1)，态密度 ν₂(E)∝E（a=4.95）→ n∝μ²（p=1.973）。纠正「缺态密度」：ν 就是 S² Dirac | 数值 |
| `exp_superconducting_doping_finiteness.py` | 超导（掺杂挂哪个「有限性」量）| 框架 λ_min=π²/N²、λ_c=2；f₀=λ_c−λ_min、f₂=ln(λ_c/λ_min) 是「谱量」（→∞）非「有限性」（→0）；1/λ_c 才是有限性；量子化 N∝λ_min^{-1/2} 是 N↔λ_min 非 n↔λ_c | 符号+数值 |
| `exp_superconducting_doping_lambda_mod.py` | 超导（N↔λ_c 两面 = N↔λ_mod）| λ_mod=log(N²/(π²ln(2N²/π²)))≈5.32，e^{λ_mod}=204.75≈207=m_μ/m_e；λ_mod≠λ_c=2（差 2.66×）——「观察者有限性」连续面是 λ_mod 非 λ_c=2 | 数值 |
| `exp_superconducting_four_faces.py` | 超导（四个面是否一个东西的投影）| λ_min、λ_mod 随 N 变（是 N 的投影），λ_c=2 不随 N 变（独立常数=经典极限）——「观察者有限性」= N（量子化）+ λ_c=2（经典极限）两独立部分 | 数值 |
| `exp_superconducting_dome.py` | 超导（穹顶/固定 Tc 矛盾）| 「掺杂=λ_c/λ_min」=3320 大数（非填充比例 0~1）；N=128 框架固定常数（无「N 随掺杂变」机制）；「掺杂=观察者有限性」且 λ_c=2 固定 → Tc 固定（无穹顶），对不上铜氧化物穹顶 35K→133K→降 | 数值 |
| `exp_topological_charge_discrete.py` | 引力侧（拓扑荷 = 离散整数 = 数）| 绕数 f=e^{ink} → n∈ℤ 离散；Hopf 荷复合律 H(g∘f)=(deg g)²H(f) 离散乘法；「粒子数=拓扑荷（离散）」框架内成立（Hopf 荷=谱流，Exp6a）——用户「物质=拓扑缺陷→数」前半坐实 | 数值 |
| `exp_modular_flow_metric.py` | 引力侧（自指不动点第②问：模流→度规）| 模流（尺度不变 φ=-log r）→ 共形平坦 R=0（∇²(-log r)=0 数值 1.5e-4）；「模流→度规」给平（退回 e=δ）不给弯曲——确认「长程 ⟂ 弯曲」墙（长程=平，弯曲=短程），第②问不破墙 | 数值 |
| `exp_vortex_long_range.py` | 引力侧（拓扑涡旋的长程性）| 2D 涡旋相互作用 V(r)=ln r（对数长程，非 1/r）；涡旋是标量（绕数∈ℤ）、度规是张量（g_μν）——「标量涡旋 → 张量度规」是「长程 ⟂ 弯曲」墙的真正形状 | 数值 |
| `exp_lambda_c_curvature.py` | 引力侧（λ_c 位置依赖→曲率）| 位置依赖 λ_c（高斯鼓包）→ 短程弯曲（R≠0 峰值 0.541，远处衰减 1.5%）=「弯曲（短程）」半（缺陷=角亏，Q1 已解），不破「长程 ⟂ 弯曲」 | 数值 |
| `exp_power_law_curvature.py` | 引力侧（幂律梯度→长程标量弯曲）| 尺度协变 Ω=r^α → R=-(d-1)(d-2)α(α+2)r^{-2α-2}（幂律长程，符号+数值 p≈3.96）；纠正「长程=平」只对尺度不变成立，尺度协变给长程标量弯曲（spin-0，Weyl C=0）。⚠️ 2026-09-30 纠错：|∇φ|² 项系数原写 (d-2)、正确是 (d-2)/2（Wald 附录 D 共形变换公式），前置系数由 -2(d-1)(d-2)α(α+1) 修正为 -(d-1)(d-2)α(α+2)，幂次 r^{-2α-2} 不变 | 符号+数值 |
| `exp_scalar_vs_tensor.py` | 引力侧（标量 vs 张量）| 共形 g=Ω²η → Weyl=0（标量=温度）；各向异性 g=diag(-Ω²,1,1,1) → Weyl≠0（张量=风向）。「标量→张量」=引入方向（vielbein e_μ^a 的 a 指标，从内部 S²/Hopf 来）| 符号 |
| `exp_q_position_dependent.py` | 引力侧（张量层照抄标量层破法）| 各向异性 g=diag(-1,1,1,a(x)) 位置依赖 → Weyl≠0（张量弯曲）；方法论：标量 Ω(r)∝r^α→R≠0、张量 q(r)→Weyl≠0，同一个「常数→位置依赖幂律」| 符号+数值 |
| `exp_qdeform_position_weyl.py` | 引力侧（q 位置依赖→Weyl≠0，坐实）| 框架真实 q 变形度规 G_ab(q)=tr_q(X_aX_b)（各向异性 G_zz≠G_xx + 扭转 G_xy）做成 q(r) 位置依赖 → 各向异性/扭转位置依赖 → 非共形平坦 → Weyl≠0（张量）。破「长程 ⟂ 弯曲」自旋 2 层（概念上）| 符号+数值 |
| `exp_weyl_long_range.py` | 引力侧（Weyl 长程性坐实）| 1-G_zz/G_xx ∝ k^{-p}，p=1.93≈2（幂律，= 量子维度 [2]_q=2cos(π/(k+2)) 的 O(1/k²) 修正）；k(r)∝r^β → Weyl ∝ r^{-βp-2}（幂律长程）——「长程+张量」共存，自旋 2 层真破（不只概念，长程性坐实）| 数值 |
| `exp_topological_charge_density.py` | 引力侧（拓扑荷密度钥匙合法）| 拓扑荷（绕数）离散 ∈ℤ、拓扑荷密度（绕数/面积）连续（实数）——「物质=拓扑荷密度」合法（物质密度需连续）；但「k(r)=拓扑荷密度」仍是物质源问题（观察者内生物质，猜测非推导）| 数值 |
| `exp_three_bridges.py` | 引力侧（三座桥 = 标量线 + 张量线）| 桥三 ν₂→n（标量，n∝μ² 已通）、桥一 ν₂→k（标量，照抄桥三，缺「敏感度」权重）、桥二 ν₂→T_μν（张量，不能照抄——「权重=E」只给 T_00 标量，需「标量→张量」= q 位置依赖→Weyl）| 数值 |
| `exp_attention_unified.py` | 引力侧（注意力统一 w(E)=ρ）| ρ=C/λ（权重）、模流（log ρ 生成）、E 的秩（D 维数）三者由 ρ（迹 τ）统一 = 同一观察者的三面；「敏感度」权重 w(E)=ρ（不是新对象）；但 ρ=无偏好，掺杂=偏离无偏好（有偏好）| 符号 |
| `exp_fracture_strength.py` | 引力侧（路二：无偏好→断裂 = 结构选择门）| 均匀 D（无偏好）Tr(D⁴)=640、π-flux D（断裂）Tr(D⁴)=256（全局最优）——「无偏好→断裂」= 结构选择门（自发对称性破缺，框架已有）；断裂强度=π-flux 偏差（标量）；剩「断裂强度→物质密度」标量映射（物质=缺陷最后一环）| 数值 |
| `exp_fracture_to_matter.py` | 引力侧（断裂→物质 = 三步）| 断裂（π-flux）→结构（Dirac 点，固定）、掺杂（μ）→密度（n∝μ²，可变）——「断裂强度→物质密度」非直接映射，剩「断裂→掺杂」（π-flux→μ，μ 是输入）| 数值 |
| `exp_coarse_grain_fracture.py` | 引力侧（二元断裂→连续掺杂的粗粒化）| 局域二元 π-flux → 粗粒化 → 连续密度（0~1，标准统计力学）；但「连续密度=μ」+「局域断裂来源（涨落）」两步没推 | 数值 |
| `exp_final_wall.py` | 引力侧（最终墙：断裂 vs 掺杂是不同维度）| 断裂密度（无量纲，结构 0~1）vs 掺杂 μ（能量量纲，填充）——「断裂→掺杂」是「结构→填充」跨维度映射，需「涨落→断裂→粗粒化→μ」动力学（真边界，非新数学）| 数值 |
| `verify_hellmann_feynman.py` | 几何→物质桥（第 1 步：键序 = δE/δD）| Hellmann–Feynman：$K_{ij}=\delta E/\delta D_{ij}=\rho_{ij}$（密度矩阵元），数值有限差分误差 3.9e-7 | 数值 |
| `verify_bond_order.py` | 几何→物质桥（键序 rank-2）| 密度矩阵 $\rho_{ij}=\langle c_i^\dagger c_j\rangle$ 是 rank-2（对角=占据、非对角=相干），「标量→张量」假墙澄清 | 数值 |
| `derive_bond_Tmunu.py` | 几何→物质桥（第 2 步：键序→T_μν）| $T_{00}=-t\sum\rho_{ij}$（对角）、$T_{0i}=-it\sum(\rho-\rho^\dagger)$（非对角相位=电流 0.1488）| 数值 |
| `derive_Tmunu_noether.py` | 几何→物质桥（T_μν 严格推导）| Dirac 作用量 + Noether → $T^{\mu\nu}=i\bar\psi\gamma^\mu\partial^\nu\psi$ → 非相对论极限 → 键序公式；守恒 $\partial_\mu T^{\mu\nu}=0$（两项抵消）| 符号 |
| `derive_coarsegrain.py` | 几何→物质桥（第 3 步：粗粒化）| 窗口平均 + 拓扑标度律（芯<<σ<<间距）+ 拓扑荷守恒（绕数 8→8）| 数值 |
| `close_einstein.py` | 几何→物质桥（第 6 步：标量 Einstein）| 角亏 δ=4πGm → $R=2\delta=8\pi G T_{00}$（m=0.5/1/2 闭合）| 数值 |
| `close_tensor_einstein.py` | 几何→物质桥（第 6 步：张量 Einstein）| $G_{\mu\nu}=8\pi G T_{\mu\nu}$ 三分量（能量→牛顿势、动量→引力磁、应力→引力电）| 数值 |
| `derive_nonlinear_einstein.py` | 几何→物质桥（非线性 EH + 反向映射）| FRW 完整 $G_{\mu\nu}$ 给 Friedmann 方程 $3(\dot a/a)^2=8\pi G\rho$（非线性）；弯曲时空 Dirac $(i\gamma^\mu\nabla_\mu-m)\psi=0$（反向映射）| 符号 |
| `derive_Tij.py` | 几何→物质桥（T_ij 应力）| 应力=动量流 $T_{ij}=\frac1N\sum_k(\partial\varepsilon_k/\partial k_i)k_j n_k$，$T_{xx}=T_{yy}=0.8095$（压力）| 数值 |
| `verify_su2_k.py` | 几何→物质桥（自旋联络：量子化绕开分离）| Chebyshev 截断 $\delta=2\cos(\pi/(N+1))\Rightarrow\Delta_N=0$ → SU(2)_k（level k=N-1）独立于 U(1) 磁通，分离墙量子层消解 | 符号 |
| `verify_weyl_Tmunu.py` | 几何→物质桥（Weyl ≠ T_μν）| FRW 共形平坦 Weyl=0 但 $T_{\mu\nu}\neq0$（$C_{0101}=C_{1212}=C_{2323}=C_{0123}=0$）——几何≠物质 | 符号 |
| `verify_quantum_metric.py` | 陈数-涡旋对应（量子度规 ∫g）| QWZ 模型 ∫g 随 m 平滑（0.240→0.097），C 是 0/1 阶跃——∫g 模型依赖，不 ∝C 不 ∝C² | 数值 |
| `verify_chern_tc.py` | 陈数-涡旋对应（陈数 vs 载流子密度）| p+ip 模型：C 是拓扑阶跃、n 是光滑费米海——两者独立，f(C)∝C² 不从体来 | 数值 |
| `verify_large_C_fixed.py` | 陈数-涡旋对应（大陈数 ∫g 标度）| 手征 p+ip 缠绕 N（m=-1 修复陈数 bug）：∫g/C≈0.113 常数 → ∫g ∝ C（线性），不是 C² | 数值 |
| `verify_kramers_robust.py` | 内禀 Kramers（对 T 破缺鲁棒）| π-flux T'²=−1 对通用 T 破缺（t₂）破简并（[1,1,...]）——「鲁棒」是特定反对易类，非通用磁杂质 | 数值 |
| `exp_kramers_pairing_class.py` | 涌现 Kramers → DIII 类（第1步）| 实空间 J 与键配对的对易判据（混合，非局域 J 的教训）| 数值 |
| `exp_kramers_pairing_4x4.py` | 涌现 Kramers → DIII 类 | 4×4（子格×谷）配对分类，DIII vs D 分岔（符号坐实）| 符号 |
| `exp_kramers_commute_check.py` | 涌现 Kramers → DIII 类 | J 与 s/p 配对的对易（实空间，混合）| 数值 |
| `exp_kramers_bdg_survival.py` | 涌现 Kramers → DIII 类 | T' 存活判据（实空间，全偶重数）| 数值 |
| `exp_kramers_majorana_count.py` | 涌现 Kramers → DIII 类 | 涡旋芯 Majorana 计数（谷内 vs 谷间）| 数值 |
| `exp_kramers_momentum_class.py` | 涌现 Kramers → DIII 类 | 动量空间判据（干净，钉死 s→DIII / p→D）| 数值 |
| `exp_diii_full_gap_scan.py` | 涌现 Kramers → DIII 类 | 完全 gap 扫描（Bogoliubov-Fermi 表面 cos²kx+cos²ky=0.4）| 数值 |
| `exp_diii_full_gap_construct.py` | 涌现 Kramers → DIII 类 | 完全 gap 构造（谷三重态 d 矢量）| 数值 |
| `exp_diii_valley_triplet.py` | 涌现 Kramers → DIII 类 | 谷三重态 τ_x/τ_z（{τ,J}=0，谷 SU(2)）| 符号 |
| `exp_valley_triplet_construct.py` | 涌现 Kramers → DIII 类 | 谷三重态精确构造（非局域）| 数值 |
| `exp_diii_valley_triplet_pairing.py` | 涌现 Kramers → DIII 类 | 谷三重态配对 | 数值 |
| `exp_diii_dvector_final.py` | 涌现 Kramers → DIII 类 | Δ=σ_y⊗(d·τ) 费米统计 + T' 存活（符号坐实，3He-B 对应）| 符号 |
| `exp_diii_dvector_gap.py` | 涌现 Kramers → DIII 类 | τ 反对易 P，T' 存活 rel=0.000（实空间坐实）| 数值 |
| `exp_diii_dvector_momentum.py` | 涌现 Kramers → DIII 类 | 动量奇 d(k) 在 Γ 点 d(0)=0，Bogoliubov-Fermi 表面拓扑必然 | 数值 |
| `exp_doping_DD_selfadjoint.py` | μ 墙破解（D-D 自反）| 掺杂 = D 与 D' 自反耦合 + 谱不对称 → 费米能级差 → 粒子转移（连续掺杂）| 数值 |
| `exp_doping_interface.py` | μ 墙破解（接口）| μ_eff=ΔEF → T_BKT=μ_eff/32 欠掺升 → min(T_BKT,T_mf) 出穹顶 | 数值 |
| `exp_doping_calibration.py` | μ 墙破解（量级标定）| 杂质浓度→转移量→μ_eff→T_BKT 量级 | 数值 |
| `exp_doping_calibration_exact.py` | μ 墙破解（精确标定）| 固定粒子数 μ(p)，对接真实铜氧化物卡强关联 | 数值 |
| `exp_entanglement_mu.py` | μ 墙（纠缠谱）| 纠缠谱找 μ（十八轮判负之一）| 数值 |
| `exp_entanglement_mu_scan.py` | μ 墙（纠缠谱扫描）| 纠缠谱 μ 扫描 | 数值 |
| `exp_strong_correlation.py` | 强关联重整化 | 缺陷→曲率 / 自指→自能 / 关联群三路判负，收敛付费桥 2 | 数值 |
| `exp_doping_self_energy.py` | 强关联（思路一）| D-D 自反→自能 Σ=t_p²G₂（一体杂交自能，负结果）| 数值 |
| `exp_bridge_fee.py` | 强关联（过桥费）| 局域化动能 2t（后被「泡利禁止」修正）| 数值 |
| `exp_pairing_locality.py` | 强关联（配对局域性）| 同格点配对=0（泡利），配对只在最近邻（r=1 A-B）| 数值 |
| `exp_vortex_repulsion.py` | 排斥几何来源 | 同号拓扑荷排斥（费米子 π 磁通 Z₂ 无同号，负结果）| 数值 |
| `exp_centrifugal_repulsion.py` | 排斥几何来源 | 离心排斥=角动量动能（束缚态按角动量分层，排斥 1.07）| 数值 |
| `exp_phase_stiffness_vs_vortex.py` | 排斥几何来源（判据）| 涡旋能量=πJ_s ln(L/a)（系数 3.17≈π），排斥=相位刚度两面 | 数值 |
| `exp_bkt_duality_verify.py` | 对偶推 BKT | BKT 普适跳变 J_s(T_c)/T_c=2/π=0.6366（MC 螺旋模量）| 数值 |

## 质量谱审计 + N=128 有向区分推导链（2026-10-01，22 脚本）

> 这一轮是质量谱数值层的独立审计：发现「无自由参数推出 207」不成立（三个选择点：λ_min 选 π²、步长 c=1、N=2⁷=128 事后对齐 + 群论计数错误）。八条负路径重推 N 全撞「7=2³−1 是 3 个 qubit 状态数，不是 7 个 qubit」；但「有向区分」公设的连续面（有向→cos θ）破局，给出 N=128 的候选推导链。见 vault 笔记 [[N=128的来源：从有向区分到Chebyshev（候选推导链）]] + [[质量谱审计收口总结（三层状态+四个选择点+四条路）]]。

### 审计四路（N=128 事后对齐坐实）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_mass_route3.py` | 路 3：7→2⁷ 换其他组合量 | 只有 2⁷ 命中 207（≈反解值 128.7），事后对齐指纹 |
| `exp_mass_route4.py` | 路 4：群论计数 | Z₂³ 只有 3 自由度，7=2³−1 是状态数，正确计数 8 或 7 不给 207 |
| `exp_mass_route5_scan.py` | 第 6 步：全组合扫描 | N∈{3,7,8,14,21,49,128,168,5040}×λ_min，唯一命中 207 是 N=128+π² |
| `exp_three_Z2_uniqueness.py` | 三个 Z₂ 唯一性（缺口10） | 观察者侧 s 唯一 = Aut(ℝ⁺) 连续对合自同构（群论定理） |

### 八条负路径（重推 N，不看 207）

| 脚本 | 路径 | 结论 |
| --- | --- | --- |
| `exp_neutrino_N.py` | 中微子观测反解 N | N≈142（差 19%，同量级不命中） |
| `exp_neutrino_tension.py` | 中微子张力 | 40.7 < 50.1 meV（Δm²₃₁ 下限，差 23%） |
| `exp_neutrino_tension_v2.py` | 中微子张力 v2 | Δm²₂₁ 隐藏张力 12.6 倍（代际比 207 内在矛盾，中微子=207 从真预言降级为候选+数据张力） |
| `exp_N_relative_version.py` | D-D 自反找 N | N 是结构侧，不在填充层 |
| `exp_fermi_surface_N.py` | 费米面面积 | 2D 费米面上限 4π²≈39.5 < 128 |
| `exp_spiral_area.py` | 螺旋面积 | O(1)，不给 128 幂次 |
| `exp_spiral_curvature_dof.py` | 螺旋曲率自由度导数 | 一次/二次式，不给 128 |
| `exp_doubling_steps.py` | 自指加倍 | 「7 个加倍」=「7 个自由度」同一错误 |
| `exp_selfref_recursion.py` | 自指递归/高阶自同构 | 平方递归 N→N² 给无限维，指数递归给 168 非 128 |
| `exp_power_set_128.py` | 幂集 128 | 假出口（幂集=独立 qubit=自由度） |
| `exp_power_set_constraint.py` | 幂集约束 | 7 元素不独立，合法子集 26 非 128 |

### 有向区分推导链（N=128 候选来源，正结果）

| 脚本 | 段 | 内容 |
| --- | --- | --- |
| `exp_directed_distinction.py` | 段 1-5 | 3Z₂→7灯→128模式→布尔傅里叶→有向桥 |
| `exp_real_part.py` | 段 3-4 | 有向=反厄米（±i）、取实部=自反性（对称/反对称分解） |
| `exp_half_circle.py` | 段 5-9 | cos 偶函数→半圆、Chebyshev、+1 尺度破缺 |
| `exp_fill_gaps.py` | 段 6-7 | 布尔傅里叶、无偏好→角度均匀（补齐跳步） |
| `exp_fix_modular.py` | 纠错 | 纠正「θ=log λ」误用模流（模流=时间维度，非角度），角度均匀=旋转对称 |
| `exp_circle_half_contradiction.py` | 缺口 A | 圆 vs 半圆参数化 bug（圆→65 值，半圆→128 值） |
| `exp_radix_choice.py` | 缺口 B | 数值序 vs 角度序（进制/排序） |
| `exp_close_gaps_AB.py` | 缺口 A+B 闭合 | 区分=二分递归（编号 n 单调对应角度 θ），半圆上做 |

## 29.2 结构路线探底（2026-10-02，5 脚本）

> 缺口 3「第一→第二代 29.2」：先推结构、后看数（不反推），四个候选全排除。见 vault [[族维度：轻子λ_c机制推广到夸克（第一刀：颜色三=代三S₃两面）]] 第二十六刀（26.1–26.9）。

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_gap3_29p2.py` | 数据边界：29.2 = (m_c/m_u)/(m_s/m_d) 相对误差 ~20%，27(3³)/32(2⁵) 都在 1σ 内 | 蒙特卡洛 + 一阶传播 + 不对称敏感性三重交叉 |
| `exp_gap3_signature_rotation.py` | 号差面严格性 P1–P6：转置=反射（det=−1，平方=+I）、号差 J²=−I（±i）、「转置=号差」是松散类比（非同构） | 符号 |
| `exp_gap3_transpose_factor.py` | 先推结构：转置对称化 P_sym=(I+(12))/2 迹2、反对称化 P_anti=(I−(12))/2 迹1、3-循环完全对称化迹1 → 自然因子 {3/1,3/2,2/1} 不给 29.2 | 符号 |
| `exp_gap3_power_source.py` | 幂次来源=「区分=二分」递归→2 的幂（N=128=2⁷）；3 的幂需三分递归（框架没有）→ 27=3³ 排除；三个转置生成 S₃（不独立）、转置类大小 3（状态数） | 符号 |
| `exp_gap3_s3_binary.py` | 5 混三分对象：S₃ 非平凡元素 5 = 3 转置（阶2）+ 2 个 3-循环（阶3），二分选择只作用阶2 → 32=2⁵ 排除；Z₂³ 非平凡元素 7 全阶2（对照无矛盾） | 符号 |

**四个候选全排除**：27=3³（无三分递归）/ 32=2⁵（5 混三分对象）/ 转置投影维数比（给 {3/1,3/2,2/1}）/ 缺陷能级=3 代（框架谱无 3 重，复用 `exp_generation_3_probe`）。**∴ 29.2 结构路线（表示论/二分递归/缺陷谱）原则上给不出**，需动力学机制（Yukawa 跑动/RG，v12 附录六路线之争）。**29.2 = 已探明的开放，记录位置，主线转 CKM/PMNS。**

## 质量谱动力学探底（2026-10-03，13 脚本）

> 从「29.2 具体数字」偏回「质量谱动力学」纲领，撞出几乎闭环的链。见 vault [[质量谱：时间+观察者截断+号差旋转（精确化）]] §十三 + [[四费米子相互作用：δD玻色场积掉与关联刚度软模（质量谱动力学探底）]] + [[付费桥2精确化：π磁通SU(2)是动量空间的（不可约分解）]]。

### 完整推导链 + 公式（质量谱动力学的核心）

> 一句话：**质量层级（指数）来自「穿过有限性」（隧穿），指数里的 n = 量子化整数（= 代阶 = 绕数 = 链接数 = 阶）**。全程程序核实（sympy 群论 + 独立推导 + Jones 精确计算）。

**推导链（每环程序核实）**：

$$\underbrace{\text{公设「观察 = 一次有向区分」}}_{\text{区分 = 二分}} \xrightarrow{J^2=-I} \underbrace{\text{2 值（±i）}}_{\text{有限}} \xrightarrow{\text{2 值同构}} \underbrace{\text{泡利 }\{0,1\} + \text{截断 }[\lambda_{\min},\lambda_c]}_{\text{有限性两面}} \xrightarrow{\text{指数截断}} \underbrace{\rho=\tfrac1\lambda e^{-\lambda/\lambda_c}}_{\text{隧穿}} \xrightarrow{\lambda_{\text{mod}}=\log\frac{C}{\lambda_{\min}}} \underbrace{m_n=e^{\lambda_{\text{mod}}n}}_{\text{质量指数}} \xrightarrow{n=\text{阶}} \underbrace{n=1,2,3}_{\text{量子化整数}}$$

**关键公式**：

| 环节 | 公式 | 脚本 |
| --- | --- | --- |
| 公设 | 观察 = 一次有向区分，$J=\begin{pmatrix}0&1\\-1&0\end{pmatrix}$，$J^2=-I$ | `exp_mass_once_axiom` |
| 二分→有限 | $J^2=-I \Rightarrow$ 本征值 $\pm i$（2 值）$\Rightarrow$ 有限（2 是有限整数） | 同上 |
| 2 值同构 | 区分 ±i ↔ 泡利 {0,1} ↔ 截断 {λ_min,λ_c}（同一个「一个 bit」的三个实现） | `exp_mass_once_axiom2` |
| 隧穿（观察者态） | $\rho=C/\lambda$，$C=1/\ln(\lambda_c/\lambda_{\min})$，指数截断 = 隧穿概率幅 | `exp_mass_principle_once` |
| λ_mod（模流频率） | $\lambda_{\text{mod}}=\log\frac{C}{\lambda_{\min}}=\log\frac1{\lambda_{\min}}-\log\ln\frac{\lambda_c}{\lambda_{\min}}$（IR 主导 7.41，UV 次领头 2.09，N=128 时 λ_mod=5.32） | `exp_mass_lambda_mod_structure` |
| 质量指数 | $m_n=e^{\lambda_{\text{mod}}n}$ | `exp_mass_winding_generation` |
| 代阶 | 两个 Z₂ → S₃ → 共轭类阶 {1,2,3} | `exp_mass_winding_verify` / `exp_mass_link_order` |
| 链接数=阶（经典极限） | 环数（cycle 数）= n − 阶 + 1（阶1→环3、阶2→环2、阶3→环1） | `exp_mass_jones_order` |
| 链接数=阶（量子精确） | Kauffman 括号最低 d 次幂 = n − 阶（阶1→d²、阶2→d¹、阶3→d⁰） | `exp_mass_jones_order_quantum` |

**程序核实状态**：凡「数学结构 / 结构来源 / 群论对应 / Jones 精确 / 2 值同构」全部 sympy/numpy 坐实；只剩 2 个本体论/类比（隧穿=有限性、2值为何是占据/尺度）如实标注「本质无法程序化」。

**⚠️ 未精确接回**：质量指数 $m_n=e^{\lambda_{\text{mod}}n}$ 的形式有了、n=阶坐实了，但 $e^{\lambda·n}$ 到 207/3477 的具体数值仍是候选（之前审计「无自由参数不成立」——三个选择点 λ_min/步长/N=128 未独立逼出，见 [[质量谱审计收口总结（三层状态+四个选择点+四条路）]]）。

### 29.2 探底（5 脚本）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_gap3_generation_mapping.py` | 2−δ_N/λ_mod 代分量比能否先验给 29.2：2−δ_N 是 O(1/N²) 精细小量、候选「代→N」映射跨代比全远小于 577/19.8、反解命中 N_down≈16 是事后合理化 | 数值 |
| `exp_gap3_color_triadic.py` | 颜色三分→三分递归：3-循环=两转置乘积、阶3≠分支数3（「两个不同的3」）、dim(3⊗3⊗3)=27 是维数非递归 | 符号 |
| `exp_gap3_vacuum_propagator.py` | 背景传播子 G₀ 对角=0（真空不空只在非对角）；π-flux 4 零模（2 Dirac 点×2）=掺杂自由度（断裂→掺杂结构侧） | 数值 |
| `exp_gap3_vacuum_correlation.py` | 真空关联 G₀ 随距离衰减：异子格 ~1/r^1.6 幂律长程 | 数值 |
| `exp_gap3_density_correlation.py` | 密度-密度关联 χ=-|G₀|²~1/r^3.2（排斥，泡利反关联，Wick 定理） | 数值 |

### 攻付费桥2（4 脚本）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_bridge2_local_U.py` | 缺陷局域 SU(2) 重叠→局部U：δF 精确局域 2 矩阵元（δ函数式）、两缺陷精确不相交、翻转 edge 破坏 J_p²=−I | 数值 |
| `exp_route_A_scale.py` | 路A 真空长程关联给质量层级：幂律（相互作用）≠指数（质量层级），标度正交 | 数值+符号 |
| `exp_route_B_valley.py` | 路B 两个谷作局域双通道：零模 PR≈64-80、主导动量在 Dirac 点（动量空间非局域） | 数值 |
| `exp_route_C_sw.py` | 路C 两 D 耦合+SW：给非局域 superexchange（J=4t²/Δ），不给局域 U | 解析 |

### 质量谱动力学地基 + 绕数=代（4 脚本）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_mass_tunneling_from_interaction.py` | WKB 隧穿 λ=0（问错对象：隧穿穿的是「公设的有限性」，不是 WKB 势垒） | 数值（负） |
| `exp_mass_principle_once.py` | 公设「一次」=有限性→泡利（n∈{0,1}）+截断（λ∈[λ_min,λ_c]），三者都是二分/有限 | 符号 |
| `exp_mass_lambda_mod_structure.py` | λ_mod 独立推导：log(C/λ_min)=log(1/λ_min)−log(ln(λ_c/λ_min))，IR 主导+UV 次领头 | 符号+数值 |
| `exp_mass_winding_generation.py` | 绕数=代=量子化整数（S¹ 绕数 0,1,1 ≠ 阶 1,2,3，绕数是拓扑荷=量子化整数），闭环 | 符号 |

### 绕数=代=链接数=阶 程序核实（4 脚本）

> 把「绕数=代」从「身份识别」补到「程序核实」——群论 + Jones 精确计算。见 vault [[质量谱：时间+观察者截断+号差旋转（精确化）]] §十三。

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_mass_winding_verify.py` | 绕数=代结构来源：两个 Z₂（转置，阶2）生成 S₃（6 元素）、S₃ 共轭类阶 {1,2,3} | 符号（群论） |
| `exp_mass_link_order.py` | 链接数=阶：共轭类代表元的阶 {1,2,3}（不是类大小 {1,3,2}） | 符号（群论） |
| `exp_mass_jones_order.py` | 链接数=阶 经典极限：环数（cycle 数）= n − 阶 + 1（阶1→环3、阶2→环2、阶3→环1） | 符号（Permutation） |
| `exp_mass_jones_order_quantum.py` | Jones 量子值最低 d 次幂 = n − 阶（阶1→d²、阶2→d¹、阶3→d⁰），Kauffman 括号 d 次幂按阶分层 | 符号 |

### 公设→有限→泡利+截断 程序核实（2 脚本）

> 把「公设→有限性」从「诠释桥」坐实为「数学」——见 vault [[质量谱：时间+观察者截断+号差旋转（精确化）]] §十三。

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_mass_once_axiom.py` | 公设→有限：区分=二分（J²=−I→±i→2值→有限）是数学，不是诠释桥 | 符号 |
| `exp_mass_once_axiom2.py` | 有限→泡利+截断：2 值同构（区分 ±i、泡利 {0,1}、截断 {λ_min,λ_c} 是同一个「一个 bit」的三个实现） | 符号 |

## 缺口 8 八刀 + 缺口 11 1生2（2026-10-02，2 脚本）

> 缺口 8（电弱标度 $v\sim M_P/N_{\text{int}}^8$ 的「8」）查证 = 巧合；缺口 11（Y_ν 动力学 / 1生2）形式化精确化。

### 缺口 8 八刀查证（`exp_gap8_badao.py`）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_gap8_badao.py` | 「8」是事后对齐：$\log_{128}(M_P/v)=7.92$ 非精确 8，k=8 给 169 GeV vs 246 GeV（差 1.45 倍），「8」只是 7/9 间最近整数；三候选全排除（$2^3$ 群阶≠指数、Chebyshev N=8 是 N 值≠指数、$2\times4$ 无机制） | 数值+符号 |

**结论**：$v\sim M_P/N_{\text{int}}^8$ 的「8」= 巧合（事后对齐），非结构来源。「群阶 8 ≠ 指数 8」是对象错误（AGENTS 三问①）。$v$ 绝对标度属「标定输入」档。推荐优先级变 **11 > 10 > 12**（8 关闭）。

### 缺口 11 1生2 形式化（`exp_gap11_Ynu_dynamics.py`）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_gap11_Ynu_dynamics.py` | 精确分界：【严格】$m=\Lambda_{\text{eff}}\lambda_{\min}^{1/2}$（能隙）+ seesaw 平方 $m_D^2/M_R$（矩阵对角化）；【候选】「区分次数=$Z_2$个数 $k$→$\lambda_{\min}^{k/2}$」是相对化模式假设，非公理推导 | 符号 |

**结论**：「1生2」从「有理由的猜测」精确化为「候选链」，核心一步（区分次数→λ_min 幂次）是【未证实/待推导】档，非【否证】。与缺口 15 的 N=128 推导链是同一「候选链」模式。见 vault [[尺度读出与宇宙学预言：f 来源完整推导]] §9.8。

## 超导 Tc 线 + 色禁闭线（2026-10-02，27 脚本）

> 从「2.5 分数维」三锚点挖穿出发，推到超导 Tc 统一公式 + 色禁闭四层。见 vault [[2.5分数维线索收口：三锚点挖穿与对标失败]] + [[超导Tc的温度来源：KMS条件与模流频率（从N=128到112K）]] + [[声子=缺陷网络的集体模：BCS的ω_D与普适Tc公式的缺口]] + [[D谱编码电子不编码原子核：超导与核质量的本体论边界]] + [[色禁闭=无外部观察者的推论：框架的独特角度]]。

### 超导 Tc 线（15 脚本，`exp_sctc_*`）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_sctc_lambda_mod.py` | λ_mod 修正：Tc=Λ/64 偏大 3.3 倍 → λ_mod=5.32 偏 1.2 倍（模流频率替代 λ_c=2） | 数值 |
| `exp_sctc_lambda_mod2.py` | λ_mod 修正的脚本修正（f-string 转义 bug） | 数值 |
| `exp_sctc_lambda_sensitivity.py` | λ_mod 对 N 的敏感性：e^λ_mod=207 只在 N=128 撞，N 微调破坏 | 数值 |
| `exp_sctc_lambda_mod_exact.py` | λ_mod = log ρ 谱最大特征值，精确等式（差=0），log log 修正 = 归一化常数 C | 数值 |
| `exp_sctc_kms_temp.py` | KMS 条件：温度 T=Λ/λ_mod（模流频率 = 逆温度×能量），Tc=T/32=112K | 数值 |
| `exp_sctc_N_knob.py` | N 是温度旋钮：Tc(N) 单调降，室温 300K 对应 N≈17 | 数值 |
| `exp_sctc_phonon.py` | 声子 ω_D = J_s·λ_mod = Λ/16π = 381K，对标德拜温度 | 数值 |
| `exp_sctc_lambda_bcs.py` | λ_BCS = 2/π（BKT 跳变常数），Tc_BCS = 1.13·ω_D·e^{-π/2} = 89K | 数值 |
| `exp_sctc_2overpi.py` | 2/π 结构来源 = 涡旋能量 πJ_s ln ↔ 涡旋熵 2 ln 对偶平衡，2/π = 维度/π | 数值 |
| `exp_sctc_unified_fit.py` | 统一公式 T_c≈54Λ(eV)K 对标：铜氧化物偏 0.81~0.97，Nb 失效 | 数值 |
| `exp_sctc_lambda_dual.py` | 三问①：质量(指数 e^λ_mod) vs 温度(倒数 1/λ_mod) 是观察者截断的两个反相关投影 | 数值 |
| `exp_sctc_tunneling_thermal.py` | 三问①到底：指数 vs 倒数 = 隧穿(量子)↔热激活(统计)对偶，同源(观察者势垒 λ_mod) | 数值 |
| `exp_sctc_hydride_fit.py` | 补 M 输入：氢化物 Tc 对标，λ≈0.47(中等耦合)非 2/π(强耦合)，两条室温路 | 数值 |
| `exp_sctc_webfit.py` | 对标实测：强耦合铜氧化物偏 1.02~1.22 / 中等耦合氢化物偏 0.94~0.99 / 1层铜氧化物+弱耦合失效 | 数值 |
| `exp_sctc_order_parameter.py` | 序参量对偶 → λ_BCS=2/π 结构同源（幅度=BCS/相位=BKT，对偶平衡 2/π=维度/π） | 数值 |

### 色禁闭线（7 脚本，`exp_confinement_*`）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_confinement_su3k.py` | SU(3)_k 量子维度 [3]_q = 3-4sin²(π/(k+3))，k=2 给 φ（黄金比） | 数值 |
| `exp_confinement_mp_ratio.py` | 核质量比 m_p/m_e 两候选：N²/π²=1660（差10%）/ 6π⁵=1836（差0.001%） | 数值 |
| `exp_confinement_lambda_qcd.py` | Λ_QCD = 尺度破缺×电弱标度 = 148 MeV（反解凑），卡色的 β 函数 | 数值 |
| `exp_confinement_seesaw.py` | Λ_QCD 的 seesaw 式幂次 8 候选（Λ_eff×λ_mod=110 / (m_e·m_μ·v)^(1/3)=237 等），全反解凑没命中 200 MeV | 数值 |
| `exp_confinement_singlet.py` | 结构禁闭走通：3⊗3⊗3 = 10⊕8⊕8⊕1（色单态 = 重子 = 三夸克），群论 | 数值 |
| `exp_confinement_flip.py` | 翻转：Λ_QCD 绝对数值 → 结构因子（Λ_QCD/m_e=391 等），结构因子没找到干净组合 | 数值 |
| `exp_confinement_goldstone.py` | GMOR 锚点：手征破缺→Goldstone 介子 m_π²=(m_u+m_d)⟨ψψ⟩/f_π²=12922→m_π=114 MeV，绝对卡手征凝聚 | 数值 |

### 2.5 分数维（5 脚本，`exp_2p5_*`，负结果）

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_2p5_theorem4.py` | 定理 4 E∝N^{2.5} 是平均场近似，精确 N^{1+1/d} ln N（斜率 1.75→1.5） | 数值 |
| `exp_2p5_ds25.py` | 分数谱维数 d_s → 态密度 ν∝E^{d_s-1} → n∝μ^{d_s}（d_s=2.5 链条） | 数值 |
| `exp_2p5_dispersion.py` | 色散关系纠错：抛物线 ν∝E^{d/2-1} vs Dirac ν∝E^{d-1}（Tc 标度差一个维度） | 数值 |
| `exp_2p5_layered_ds.py` | 层状格点谱维数 d_s 在 2~3 可调（L=4 给 2.55） | 数值 |
| `exp_2p5_layered_ds2.py` | 层状格点谱维数（正确 t 范围重跑） | 数值 |

### π 磁通配对（13 脚本，`exp_pi_flux_*` + `exp_bdg_*`，接口非独有，2026-09-24 迁移）

> 从 `pi_flux_pairing/` 迁移（2026-10-02）。推「π 磁通 → 配对对称性」，结论「**接口非独有**」（超导是标准凝聚态，非独有预言）。见 vault [[π磁通配对推导：谷奇配对与手征p+ip拓扑超导（接口非独有）]] + 预印本 1.17。

| 脚本 | 内容 | 验证 |
| --- | --- | --- |
| `exp_pi_flux_dirac.py` | π 磁通 Bloch 哈密顿量 Dirac 谱 + 手征 Γ + 涌现时间反演 | 符号 |
| `exp_bdg_pairing_parity.py` | BdG + 粒子-空穴对称 ⟹ 配对宇称约束 + 子格间 iσ_y | 符号 |
| `exp_pi_flux_lattice.py` | 全格点 π 磁通 D：Dirac 零模 + T²=−1 | 数值 |
| `exp_pi_flux_dos.py` | π 磁通 DOS ∝ \|E\| + 掺杂 N(μ)∝\|μ\| | 数值 |
| `exp_pi_flux_tc.py` | 掺杂 Dirac 的 Tc：BCS 型 + 量子临界阈值 V_c | 数值 |
| `exp_pi_flux_vhs.py` | van Hove 奇点：DOS 对数发散 + Tc 增强 | 数值 |
| `exp_pi_flux_chern.py` | 手征 p+ip 的 Chern 数 C=1 | 数值 |
| `exp_pi_flux_majorana.py` | 手征 p+ip 涡旋芯 Majorana 零模 | 数值 |
| `exp_pi_flux_bkt.py` | 超流刚度 J_s + BKT 温度 | 数值 |
| `exp_pi_flux_tc_3d.py` | 三维 BCS Tc（V_eff 是自由参数） | 数值 |
| `exp_pi_flux_spontaneous.py` | 自发 vs 手放 π 磁通：完全磁通平带是小图（N≤16）伪影 | 数值 |
| `exp_pi_flux_temperature.py` | 温度能否区分自发 vs 手放：ΔF/N→0 消失 | 数值 |
| `exp_pi_flux_chern_boundary.py` | 拓扑荷（陈数）能否救：体陈数对 holonomy 不敏感 | 数值 |

**结论（接口非独有，勿回潮）**：超导线推到底 = 标准凝聚态复述（配对对称性=费米统计、拓扑=Read-Green/Sato、Tc 推不出，V_eff 是接入）。根因：理论推「几何+联络+对称」，推不出「物质相互作用」。

## 运行

```powershell
py -m experiments.exp_step1_ontology_R
```

每个脚本跑完写 `experiments/<name>_last_run.json`（数值快照）。

## 完整脚本清单（242 个 = 202 + 27 超导 Tc/色禁闭 + 13 π 磁通配对）

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
- **缺口系列（18）**：`exp_gap1_seesaw` + `exp_gap1_Ynu_derive` + `exp_gap1_Ynu_path3`（缺口 1 seesaw，Y_ν=1）、`exp_gap5_wz` ×3（缺口 5 w(z) 双重负结果）、`exp_gap6_pi_factor` + `exp_gap6_pi_v2`（缺口 6 π 因子面积单位约定）、`exp_gap7_two_N`（缺口 7 两个 N 精确关系）、`exp_gap8_badao`（缺口 8 八刀查证=巧合）、`exp_gap11_Ynu_dynamics`（缺口 11 1生2 形式化精确化）、`exp_gap15_observer_independence`（缺口 15 HS 正交 7 vs 独立本征值 3）、`exp_gap15_result_space`（缺口 15 结果空间 2³）、`exp_gap15_observer_state`（缺口 15 行为空间 2⁷，对易≠可分离）、`exp_gap15_bisection`（缺口 15 区分=二分=谱分解）、`exp_gap15_lambda_min`（缺口 15 λ_min=尺度破缺 vs 边缘间隔）、`exp_gap15_step_c1`（缺口 15 步长 c=1 = 区分二分 ±1 的单位，非无绝对尺度循环）、`exp_gap15_unified_convention`（缺口 15 三个选择点合并成一个小约定：一次区分=独立操作，同时给 N=128/c=1/λ_min）、`exp_gap11_N128_convention`（缺口 11 与 N=128 共享「区分可分离」小约定、但「区分→幂次」规律不同：2^k vs λ_min^{k/2}）、`exp_gap11_distinction_discretization`（缺口 11 区分=离散化三步链条：区分=离散化→尺度破缺→质量，②定义+③④⑤严格）（见上方表）。
- **质量谱缺口 3 + 边界验证（4）**：`exp_transpose_reflection`（结构预言：转置不投影）、`exp_first_second_generation`（缺口 3 29.2 探索）、`exp_crosscheck`（预印本 1.15/1.16 交叉核对）、`exp_verify_prediction`（结论边界验证，防过度声称）（见上方表）。
- **桥 C（2）**：`exp_bridgeC_Gz`（G(z) 跑动重新检查：符号反转但量级差 62 倍 + Λ 双重角色）、`exp_bridgeC_Lambda_dual`（Λ 双重角色：λ_c vs λ_min 出路不成立，真正出路需 Λ_G≠Λ_ρ 新结构）（见上方表）。
- **超导（11，范式反转：相互作用=关联）**：`exp_superconducting_tc_endogenous`（T_BKT=μ/32 只依赖相位刚度，不依赖 V_eff + 两个刚度统一）、`exp_superconducting_tc_mu_lambdac`（掺杂 μ↔观察者截断 λ_c 候选）、`exp_superconducting_tc_room_temp`（室温无上限，Λ≈1.65eV 材料工程目标）、`exp_superconducting_tc_derivation_symbolic`（J_s 推导符号验证 + μ² 坑）、`exp_superconducting_tc_stiffness_numeric`（J_s 数值验证 + 载流子密度，抓出因子错）、`exp_superconducting_dos_spectral_dim`（态密度从谱维数导出：d_s=2→ν(E)∝E→n∝μ²）、`exp_superconducting_s2_dirac_dos`（「秤」=S² Dirac）、`exp_superconducting_doping_finiteness`（掺杂挂哪个「有限性」量）、`exp_superconducting_doping_lambda_mod`（N↔λ_c 两面 = N↔λ_mod）、`exp_superconducting_four_faces`（四个面验证）、`exp_superconducting_dome`（穹顶/固定 Tc 矛盾）（见上方表）。
- **物质侧墙（16，物质=缺陷=拓扑 → 断裂→掺杂）**：`exp_topological_charge_discrete`（拓扑荷离散=数）、`exp_topological_charge_density`（拓扑荷密度连续）、`exp_modular_flow_metric`（模流→平）、`exp_vortex_long_range`（涡旋 ln r 长程）、`exp_lambda_c_curvature`（λ_c 鼓包→短程弯曲）、`exp_power_law_curvature`（尺度协变→长程标量弯曲）、`exp_scalar_vs_tensor`（标量 vs 张量）、`exp_q_position_dependent`（张量层照抄标量层）、`exp_qdeform_position_weyl`（q 位置依赖→Weyl≠0）、`exp_weyl_long_range`（Weyl 长程性 p≈2）、`exp_three_bridges`（三座桥 ν₂ 源）、`exp_attention_unified`（注意力统一 w(E)=ρ）、`exp_fracture_strength`（断裂强度 Tr(D⁴) 640→256）、`exp_fracture_to_matter`（断裂→物质三步）、`exp_coarse_grain_fracture`（二元→连续粗粒化）、`exp_final_wall`（最终墙：断裂 vs 掺杂跨维度）（见上方表）。
- **几何→物质桥（15，键序→T_μν→Einstein→曲率，六步闭环）**：`verify_hellmann_feynman`（键序=δE/δD，Hellmann-Feynman）、`verify_bond_order`（键序 rank-2，「标量→张量」假墙）、`derive_bond_Tmunu`（键序→T_μν）、`derive_Tmunu_noether`（T_μν Noether 严格推导 + 守恒律）、`derive_coarsegrain`（粗粒化 + 拓扑标度律 + 拓扑守恒）、`close_einstein`（标量 Einstein）、`close_tensor_einstein`（张量 Einstein 三分量）、`derive_nonlinear_einstein`（非线性 EH/Friedmann + 反向映射）、`derive_Tij`（T_ij 应力=动量流）、`verify_su2_k`（SU(2)_k Chebyshev 截断，量子化绕开分离墙）、`verify_weyl_Tmunu`（Weyl≠T_μν，几何≠物质）、`verify_quantum_metric`（量子度规 ∫g）、`verify_chern_tc`（陈数 vs 载流子密度）、`verify_large_C_fixed`（大陈数 ∫g∝C 线性）、`verify_kramers_robust`（Kramers 鲁棒是特定反对易类）（见上方表）。
- **尺度读出中间（4）**：`exp_f_observer`（f=E 谱形式）、`exp_lambda_min`（λ_min/λ_c 两面）、`exp_selfref_two_Z2`（含源 EH 分类）、`exp_N_lambdac`（N↔λ_c 两面）（见上方表）。
- **质量谱中间（14）**：`exp_mass_frequency`、`exp_mass_spectrum`、`exp_signature_rotation`、`exp_lambda_c_exponentiation`、`exp_three_Z2_fano`、`exp_mass_final`、`exp_third_Z2_independence`、`exp_third_Z2_strict`、`exp_mass_form`、`exp_mass_tunneling`、`exp_mass_D2_tunneling`、`exp_mass_discrete_jump`、`exp_mass_step_one`、`exp_mass_prediction`（见上方表）。
- **超导线 DIII 类（14，涌现 Kramers → DIII 类配对）**：`exp_kramers_pairing_class`、`exp_kramers_pairing_4x4`、`exp_kramers_commute_check`、`exp_kramers_bdg_survival`、`exp_kramers_majorana_count`、`exp_kramers_momentum_class`、`exp_diii_full_gap_scan`、`exp_diii_full_gap_construct`、`exp_diii_valley_triplet`、`exp_valley_triplet_construct`、`exp_diii_valley_triplet_pairing`、`exp_diii_dvector_final`、`exp_diii_dvector_gap`、`exp_diii_dvector_momentum`（见上方表）。
- **μ 墙破解（6，D-D 自反掺杂）**：`exp_doping_DD_selfadjoint`、`exp_doping_interface`、`exp_doping_calibration`、`exp_doping_calibration_exact`、`exp_entanglement_mu`、`exp_entanglement_mu_scan`（见上方表）。
- **强关联 + 排斥几何来源 + 对偶推 BKT（8）**：`exp_strong_correlation`、`exp_doping_self_energy`、`exp_bridge_fee`、`exp_pairing_locality`、`exp_vortex_repulsion`、`exp_centrifugal_repulsion`、`exp_phase_stiffness_vs_vortex`、`exp_bkt_duality_verify`（见上方表）。

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
- [x] **超导范式反转 + 物质侧墙（2026-09-29）**：超导接口坐实（$J_s/T_{BKT}/n$ 三步验证，抓出 3 处手推因子错）；「掺杂=观察者有限性」给固定 Tc（卡点，留活口）。物质侧墙推到底（物质=缺陷=拓扑、长程⟂弯曲两层破法、三座桥、注意力统一），最终墙 = 断裂（二元）→掺杂（连续）跨维度映射（开放，需新动力学）。见 vault 预印本 1.17 / 1.18。
- [x] **几何→物质桥完整闭环（2026-09-30）**：键序 $K_{ij}=\delta E/\delta D_{ij}=\rho_{ij}$（Hellmann-Feynman）→ 完整 $T_{\mu\nu}$（Noether 推导 + 守恒律）→ 粗粒化 → 角亏（几何）→ Einstein（标量 + 张量，六步闭环）。关键澄清：「标量→张量」是假墙（键序 rank-2 直接完整）、Weyl≠T_μν（几何≠物质）、桥=键序不是 Weyl、f(C)∝C² 假假设、自旋联络分离墙量子层消解（SU(2)_k）。见 vault 笔记 [[几何→物质桥完整推导：键序→T_μν→Einstein→曲率（六步闭环+程序验证）]]。
- [x] **纠错：共形曲率公式系数（2026-09-30）**：`exp_power_law_curvature.py` 的共形标量曲率公式 $R=-2(d-1)e^{-2\phi}[\nabla^2\phi+(d-2)|\nabla\phi|^2]$ 第二项 $|\nabla\phi|^2$ 系数错（多 2 倍），正确为 $(d-2)/2$（Wald 附录 D 标准共形变换公式）。独立 sympy 验证（直接算度规 $g=r^{2a}\delta$ 的 Ricci 标量，d=3）：真实值 $-2a(a+2)$ 与正确公式 $-(d-1)(d-2)a(a+2)$ 一致，与旧公式 $-2(d-1)(d-2)a(a+1)$ 不符。已同步修正预印本 1.18 §2.2/§9.4 + 本 README + 脚本。**定性结论（长程幂律 $r^{-2\alpha-2}$）不变，只前置系数错**。教训：符号验证只保证「笔记与代码自洽」，不保证「公式本身对」——必须独立从度规算 Ricci 标量这类「独立验证」才能抓出 self-consistent but wrong。
- [x] **质量谱审计 + N=128 有向区分推导链（2026-10-01）**：审计坐实「无自由参数推出 207」不成立（三个选择点 + N=2⁷ 群论计数错误）；八条负路径重推 N 全撞「7 是 3 个 qubit 状态数」；「有向区分」公设的连续面破局，给出 N=128 候选推导链（10 段全程序验证，无新公设）。见 vault [[N=128的来源：从有向区分到Chebyshev（候选推导链）]] + [[质量谱审计收口总结（三层状态+四个选择点+四条路）]]。诚实定位：候选推导链（每步数学验证），非严格证明。
- [x] **29.2 结构路线探底（2026-10-02）**：缺口 3「第一→第二代 29.2」先推结构后看数，四个候选全符号排除（27=3³ 无三分递归 / 32=2⁵ 混三分对象 / 转置投影维数比 {3/1,3/2,2/1} / 缺陷能级=3 代谱无 3 重）。29.2 结构路线原则上给不出，需动力学机制（Yukawa 跑动/RG），主线转 CKM/PMNS。5 脚本：`exp_gap3_29p2` / `exp_gap3_signature_rotation` / `exp_gap3_transpose_factor` / `exp_gap3_power_source` / `exp_gap3_s3_binary`（+ 复用 `exp_generation_3_probe`）。见 vault [[族维度：轻子λ_c机制推广到夸克（第一刀：颜色三=代三S₃两面）]] 第二十六刀。
- [x] **缺口 8 八刀查证 + 缺口 11 1生2 形式化（2026-10-02）**：缺口 8「$v\sim M_P/N_{\text{int}}^8$ 的 8」查证 = 巧合（事后对齐，$\log_{128}=7.92$ 非 8，三候选全排除，「群阶 8≠指数 8」对象错误）；缺口 11「1生2」从「有理由的猜测」精确化为「候选链」（【严格】能隙+seesaw 平方；【候选】区分次数→λ_min 幂次，非公理推导）。2 脚本：`exp_gap8_badao` / `exp_gap11_Ynu_dynamics`。
- [x] **超导 Tc 线 + 色禁闭线（2026-10-02）**：从「2.5 分数维」三锚点挖穿出发，推到超导 Tc 统一公式 + 色禁闭四层。**① 超导 Tc**：λ_mod 修正（Tc=Λ/64 偏大 3.3 倍 → λ_mod=5.32 偏 1.2 倍）+ KMS 条件严格化（温度 T=Λ/λ_mod，λ_mod=log ρ 谱最大特征值精确、差=0）+ 声子 ω_D=Λ/16π=381K + 用户洞察「对偶 1 分为 2」（区分可分离 = N=128 小约定 → 独立声子）→ **λ_BCS=2/π**（BKT 跳变常数 = BCS 耦合常数，结构来源 = 涡旋能量 πJ_s ln ↔ 涡旋熵 2 ln 对偶平衡，2/π = 维度/π）。统一公式 $T_c\approx54\,\Lambda(\text{eV})$ K，对标铜氧化物（强耦合 2-3 层）偏 0.81~0.97 倍，弱耦合（Nb，λ≠2/π）失效 60 倍。**② 色禁闭**：无外部观察者 → 无孤立夸克（物理直觉第一性）+ SU(3)_k 量子维度 $[3]_q=3-4\sin^2\frac{\pi}{k+3}$（k=2 给 φ）+ D-D 自反给有量纲（Λ_QCD = 闭合↔开放谱差）。**卡点**：Λ_QCD 绝对能标 = 色的 β 函数（跑动耦合 → 低能发散），框架 β 函数只做过引力 G 没做过色荷 g_color。见 vault [[超导Tc的温度来源：KMS条件与模流频率（从N=128到112K）]] + [[声子=缺陷网络的集体模：BCS的ω_D与普适Tc公式的缺口]] + [[D谱编码电子不编码原子核：超导与核质量的本体论边界]] + [[色禁闭=无外部观察者的推论：框架的独特角度]]。脚本 27 个在 `experiments/`（`exp_sctc_*` 15 个 + `exp_confinement_*` 7 个 + `exp_2p5_*` 5 个），见下方「超导 Tc 线 + 色禁闭线（27 脚本）」清单。
- [x] **π 磁通配对迁移（2026-10-02）**：原 `pi_flux_pairing/` 目录的 13 个脚本（`exp_pi_flux_*` 12 个 + `exp_bdg_pairing_parity` 1 个）迁入本库 `experiments/`，原目录删除。推「π 磁通 → 配对对称性」，结论「**接口非独有**」（超导是标准凝聚态，非独有预言，见 [[π磁通配对推导：谷奇配对与手征p+ip拓扑超导（接口非独有）]] + 预印本 1.17）。超导线的所有脚本（π 磁通配对 + 超导 Tc + 色禁闭 + 2.5 分数维 = 40 个）从此统一在一个库。
- [x] **质量谱动力学探底 + 口径方案 B + λ_min 收回 + 预印本 1.3 勘误（2026-10-03）**：探「质量谱缺的动力学」能否内生，完整推导见 vault [[四费米子相互作用：δD玻色场积掉与关联刚度软模（质量谱动力学探底）]]。净产出五条：① 口径定方案 B（模流 vs 有向区分 = 时间两个面：模流给时间维度、有向区分给号差方向，已改对齐页 #186 + 阴阳图景注记 + 本 README）；② λ_min 选择点①收回（`exp_lambda_min_precision.py`：λ_min=π²/N² 是 2−δ_N 领头阶，2−δ_N=2−2cos(π/(N+1))≈π²/(N+1)² 差 O(1/N³)，结构来源=三合一「质量=尺度破缺=2−δ_N」）；③ Yukawa 耦合内生（$S_F=\langle\psi|D|\psi\rangle\to S_{\text{int}}=\sum\delta D\,\psi^\dagger\psi$，=键序 Hellmann-Feynman）；④ 四费米子=积掉 δD（$S_{\text{eff}}=\rho K^{-1}\rho$，软模相位型=关联刚度=排斥吸引对偶统一，接超导 Tc 线；`exp_interaction_vertex.py` + `exp_interaction_vertex_true_min.py`）；⑤ 预印本 1.3 勘误（「度规有质量 ≈2e-2」撤回——D* 是鞍点 8 负本征值，+2e-2 是最低正本征值非度规质量；核心「短程接触引力」不受影响）。**诚实边界**：四费米子（费米子-费米子）≠ Yukawa（费米子-希格斯），没触及质量谱缺的 29.2（Yukawa 随代跑动，补十已判结构路线给不出）。3 脚本均自包含（内联原 gauge_emergence 的 `toroidal_D`/`basis_matrices`/`simple_hessian`）。
- [x] **质量谱动力学纲领闭环（2026-10-03 续）**：从「29.2 具体数字」偏回「质量谱动力学」纲领，撞出几乎闭环的链。① 29.2 四条结构/动力学路线全判死（`exp_gap3_generation_mapping` 2−δ_N 代分量 / `exp_gap3_color_triadic` 颜色三分 / `exp_gap3_vacuum_propagator` 背景传播子 G₀ 对角=0 + π-flux 零模=掺杂结构侧 / `exp_gap3_vacuum_correlation` + `exp_gap3_density_correlation` 真空关联 G₀~1/r^1.6、χ=-|G₀|²~1/r^3.2 排斥）；② 攻付费桥2 三路判负（`exp_bridge2_local_U` 缺陷局域 SU(2) δ函数式不重叠 / `exp_route_A_scale` / `exp_route_B_valley` / `exp_route_C_sw`）——「框架给非局域，给不了局域 U」= 付费桥2；③ 质量谱动力学地基（`exp_mass_principle_once`）：公设「一次」=有限性→泡利+截断；④ λ_mod 独立推导（`exp_mass_lambda_mod_structure`）：IR 主导（7.41）+UV 次领头（2.09）→λ_mod=5.32；⑤ 绕数=代=量子化整数（`exp_mass_winding_generation`）闭环：公设「一次」→泡利+截断→隧穿→λ_mod→质量指数 e^{λ·n}（n=绕数=代阶）。**✅ 最后一环「两个 Z₂ 链接数 = S₃ 阶 1,2,3」已坐实**（群论两个 Z₂→S₃→阶1,2,3 + 经典极限环数=n−阶+1 + Jones 量子最低 d 次幂=n−阶，`exp_mass_winding_verify` / `exp_mass_link_order` / `exp_mass_jones_order` / `exp_mass_jones_order_quantum`）。**✅ 公设→有限→泡利+截断 也坐实**（区分=二分 J²=−I→±i→2值→有限是数学，`exp_mass_once_axiom`；有限→泡利+截断是 2 值同构，`exp_mass_once_axiom2`）。**仍标本体论/类比**（本质无法程序化）：隧穿=有限性、2值为何是占据/尺度、排斥吸引对偶。见 vault [[质量谱：时间+观察者截断+号差旋转（精确化）]] §十三。

## 严格证明分层（2026-09-25 修正）

- **数值验证**（numpy）：坐实，非证明；
- **符号验证**（sympy 恒等式/级数/极限/积分/函数方程解）：对可符号化的数学**就是严格证明**；
- **不可符号化**（算子范数、谱测度、泛函分析 sup）：**引用标准定理**（文献已有）或证明助手。

rigorous 脚本：`exp_step3_state_rho_rigorous`、`exp_step4_point_rigorous`、
`exp_step5_metric_rigorous`、`exp_theorem_A_rigorous`、`exp_theorem_C_rigorous`、
`exp_theorem_D_rigorous`、`exp_reference_citations`。

**唯一真开放问题 = 定理 B 非线性精确版**（弱场→全阶），其余全部严格闭合。

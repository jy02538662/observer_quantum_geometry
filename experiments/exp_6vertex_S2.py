"""
exp_6vertex_S2.py

任务：在 6-vertex 模型上，读出「弦端点的自旋 1/2 内部自由度 → 格点 S² 方向」。

付费桥 2 八节「层 3」= 从弦端点读出每格点一个 S² 方向（弦端点 = q 变形自旋 1/2 = S²_q）。

本脚本做最核心的检查：
1. 6-vertex 的 6 种顶点配置 → 6 个「自旋方向」（用 2 个 in 箭头的连线方向）；
2. 这 6 个方向是不是「S² 的离散近似」（Bloch 球面上的点）；
3. 检查「弦」（箭头路径）和「弦端点」，读出端点的自旋方向；
4. 诚实定位：6-vertex 给「离散 S²_q 方向」还是「连续 S²」？

关键障碍要先厘清：6-vertex 的「箭头」是 in/out（Z_2 离散），
而「S² 方向」是连续（Bloch 球面）。两者是否同构是本脚本要查的。
"""

import numpy as np

# ---------------- 6 种顶点配置 → 自旋方向 ----------------
# 每个顶点 4 个边（上 u、下 d、左 l、右 r），in=+1（指向顶点），out=-1
# 约束 2 in 2 out。用「2 个 in 箭头的连线方向」定义自旋方向 n ∈ S²。
# 约定：上 in = 箭头从上方指向顶点 → 方向有 -z 分量；下 in → +z；左 in → +x；右 in → -x
def config_to_direction(u, d, l, r):
    """4 个边的 in/out（±1）→ 自旋方向 n（3 维单位向量）。

    用「2 个 in 箭头的位置连线方向」：上=(0,0,1)、下=(0,0,-1)、
    左=(-1,0,0)、右=(1,0,0)，2 个 in 位置的中点方向（归一化）。
    """
    pos = {"u": np.array([0, 0, 1]), "d": np.array([0, 0, -1]),
           "l": np.array([-1, 0, 0]), "r": np.array([1, 0, 0])}
    ins = []
    if u == 1:
        ins.append(pos["u"])
    if d == 1:
        ins.append(pos["d"])
    if l == 1:
        ins.append(pos["l"])
    if r == 1:
        ins.append(pos["r"])
    # 2 个 in 的连线方向（两 in 位置的差，若平行则取其中一个方向）
    if len(ins) != 2:
        return None
    n = ins[0] - ins[1]
    norm = np.linalg.norm(n)
    if norm < 1e-9:  # 上下 in（竖直）或左右 in（水平）：用位置本身方向
        n = ins[0]
        norm = np.linalg.norm(n)
    return n / norm


# 6 种配置（2 in 2 out）：(u, d, l, r)，in=+1, out=-1
configs = {
    "a1": (+1, -1, +1, -1),  # 上 in 下 out 左 in 右 out
    "a2": (-1, +1, -1, +1),  # 上 out 下 in 左 out 右 in
    "b1": (+1, -1, -1, +1),  # 上 in 下 out 左 out 右 in
    "b2": (-1, +1, +1, -1),  # 上 out 下 in 左 in 右 out
    "c1": (+1, +1, -1, -1),  # 上 in 下 in 左 out 右 out
    "c2": (-1, -1, +1, +1),  # 上 out 下 out 左 in 右 in
}

print("=" * 70)
print("1. 6 种顶点配置 → 自旋方向 n ∈ S²")
print()
directions = {}
for name, (u, d, l, r) in configs.items():
    n = config_to_direction(u, d, l, r)
    directions[name] = n
    print(f"  {name}: ({u:+.0f},{d:+.0f},{l:+.0f},{r:+.0f}) → n = ({n[0]:+.3f},{n[1]:+.3f},{n[2]:+.3f})")
print()

# 去重：a1 和 a2 方向相反，b1/b2 相反，c1/c2 是轴
print("=" * 70)
print("2. 6 个方向在 S² 上的分布")
print()
# 收集非零方向，检查它们是不是「S² 的离散点」
dirs = list(directions.values())
# 检查是否有重复/相反
print("  6 个方向（去重后本质有几个）:")
unique = []
for n in dirs:
    is_new = all(np.linalg.norm(n - m) > 1e-6 and np.linalg.norm(n + m) > 1e-6 for m in unique)
    if is_new:
        unique.append(n)
print(f"    本质不同的方向数 = {len(unique)}")
for n in unique:
    print(f"      n = ({n[0]:+.2f},{n[1]:+.2f},{n[2]:+.2f})")
print()
print("  这 6 个方向覆盖 S² 的哪些点？")
print("    c1 = 纯 z 轴（竖直），c2 = 纯 x 轴（水平）")
print("    a1/a2, b1/b2 = 赤道（z-x 平面）的对角方向")
print("  → 6 个方向 = S² 的【离散 6 点近似】（赤道 4 点 + z 轴 + x 轴）")
print("  → 缺 y 轴！6-vertex 的箭头只覆盖 z-x 平面，不给 y 分量")
print()

print("=" * 70)
print("3. 关键障碍：6-vertex 箭头 vs S²")
print()
print("  6-vertex 的「箭头」是 in/out（Z_2），映射到自旋方向后：")
print("    - 只有 z-x 平面（缺 y 轴）")
print("    - 6 个离散点（非连续 S²）")
print("  而「S² 方向」（Bloch 球面）是 3 维连续（含 y 轴）")
print()
print("  所以 6-vertex 的「箭头」给的是 S² 的【z-x 平面 6 点离散近似】，")
print("  不是【完整连续 S²】。")
print()

print("=" * 70)
print("4. 但弦端点（自旋 1/2 内部空间）可能有完整 S²")
print()
print("  关键区分：")
print("    - 6-vertex 的「箭头」（in/out）= Z_2，z-x 平面 6 点")
print("    - 弦端点的「自旋 1/2 内部空间」= SU(2) 自旋 1/2 = 完整 S²（含 y）")
print()
print("  付费桥 2 八节说的「弦端点 = 自旋 1/2 = S²_q」指的是后者（内部空间），")
print("  不是 6-vertex 的「箭头」（in/out）。")
print("  → 要读出 S² 方向，需要读出「弦端点携带的 SU(2) 自旋 1/2 内部自由度」，")
print("    而 6-vertex 的「箭头」本身只给 Z_2（z-x 平面），不含这个内部自由度。")
print()

print("=" * 70)
print("5. 诚实定位")
print()
print("  6-vertex 模型（箭头 in/out）直接给的是「Z_2 箭头」（z-x 平面 6 点），")
print("  不是「SU(2) 自旋 1/2 内部空间」（完整 S²，含 y 轴）。")
print("  所以「读出弦端点 S² 方向」卡在：6-vertex 的箭头没有 y 轴分量，")
print("  需要额外引入「SU(2) 内部自由度」（付费桥 2 的「自旋 1/2 弦」），")
print("  而这正是「弦 = 自旋 1/2 线」的「自旋」要额外指定的东西。")
print()

print("=" * 70)
print("诚实结论见下方分析（不写进脚本）")

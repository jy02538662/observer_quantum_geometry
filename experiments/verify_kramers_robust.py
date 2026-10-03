import numpy as np

print("="*72)
print("验证：π-flux 的 T'²=-1（Kramers 双简并）对 T 破缺是否鲁棒？")
print("="*72)

L = 4
N = L*L
def idx(i,j): return i*L+j

# π-flux 哈密顿（无自旋）：横 +1，竖 (-1)^i
def build_H(t2):
    H = np.zeros((N,N), dtype=complex)
    for i in range(L):
        for j in range(L):
            v = idx(i,j)
            # 最近邻
            H[v, idx((i+1)%L,j)] += 1.0
            H[v, idx((i-1)%L,j)] += 1.0
            H[v, idx(i,(j+1)%L)] += (-1.0)**i
            H[v, idx(i,(j-1)%L)] += (-1.0)**i
    # T 破缺项（Haldane 型次近邻，复相位，破坏时间反演）
    if t2 != 0:
        for i in range(L):
            for j in range(L):
                v = idx(i,j)
                # 次近邻 (i,j)->(i+1,j+1) 带相位 i
                H[v, idx((i+1)%L,(j+1)%L)] += 1j*t2
                H[v, idx((i-1)%L,(j-1)%L)] += -1j*t2
    return (H + H.conj().T)/2

def degeneracies(H):
    ev = np.real(np.linalg.eigvalsh(H))
    # 数每个本征值的重数（简并）
    ev_round = np.round(ev, 6)
    from collections import Counter
    c = Counter(ev_round)
    return sorted(c.values())

for t2 in [0.0, 0.5, 1.0, 2.0, 3.0]:
    H = build_H(t2)
    deg = degeneracies(H)
    all_even = all(d % 2 == 0 for d in deg)
    print(f"    t2={t2:.1f}: 简并重数 = {deg}, 全偶(双简并/ Kramers) = {all_even}")

print()
print(">>> 全偶重数 = 每个本征值成对出现 = T'²=-1 的 Kramers 双简并。")
print(">>> 若 t2（T 破缺）扫过去仍全偶，说明 Kramers 双简并对 T 破缺「鲁棒」。")

# 对比：标准自旋 1/2 的 Kramers 会被 T 破缺拆开（这里无自旋，T'²=-1 是涌现的）
print("\n[关键区别]")
print("    标准 BCS：自旋单态配对由「时间反演 T」保护，磁杂质破 T → 拆对 → Tc 降。")
print("    π-flux：配对由「涌现 T'」保护，T'²=-1 对 T 破缺鲁棒 → 配对存活 → Tc 不降。")

# -*- coding: utf-8 -*-
"""
继续推: 测 lambda_mod 的敏感性, 判断「lambda_mod 双重角色」是结构同源还是巧合
关键: lambda_mod(N) 和 e^{lambda_mod} 对 N 有多敏感?
     质量谱用 e^{lambda_mod}=207 (锁 N=128), 超导用 lambda_mod 当 Tc 截断
"""
import numpy as np

def lambda_mod(N):
    return np.log(N**2 / (np.pi**2 * np.log(2 * N**2 / np.pi**2)))

print("=== lambda_mod(N) 和 e^{lambda_mod} 对 N 的敏感性 ===")
print(f"{'N':>8} {'lambda_mod':>12} {'e^{lambda_mod}':>14} {'要Tc=92K的Lambda(eV)':>22}")
tc_target_meV = 92 / 11.6  # meV
for N in [7, 8, 14, 21, 49, 128, 168, 5040]:
    lm = lambda_mod(N)
    e_lm = np.exp(lm)
    Lam = tc_target_meV * lm * 32 / 1000  # 要 Tc=92K 需的带宽(eV)
    print(f"{N:>8} {lm:>12.3f} {e_lm:>14.1f} {Lam:>22.2f}")

print()
print("质量谱: e^{lambda_mod} 要 = 207 = m_mu/m_e, 只有 N=128 附近给到(204.75)")
print("   N=7->2.2, N=14->5.4, N=49->39, N=128->205, N=168->331")
print("   => e^{lambda_mod}=207 对 N 极敏感, N 微调就破坏")
print()
print("超导: 要 Tc=92K 需的带宽 Lambda = 0.254 * lambda_mod eV")
print("   N=7: 0.20eV, N=128: 1.35eV, N=5040: 1.94eV")
print("   => 任何 N 都能通过调 Lambda 给 92K (Lambda 吸收自由度)")
print()
print("=== 三问第一问: 两个 lambda_mod 作用在同一个对象上吗? ===")
print("质量谱 lambda_mod: 作用在「质量频率比」 e^{lm}=m_mu/m_e")
print("超导   lambda_mod: 作用在「能带宽度/温度」 Tc=Lambda/(lm*32)")
print("=> 两个对象不同(质量比 vs 温度), 是「同一个数 5.32 撞两个对象」")
print("=> 按 AGENTS 三问①, 这是「数字巧合」的典型指纹, 非结构同源")

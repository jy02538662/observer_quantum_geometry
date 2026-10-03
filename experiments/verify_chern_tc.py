import numpy as np

print("="*72)
print("数值测试：p+ip 超导里，载流子密度 n 和陈数 C 是什么关系？")
print("="*72)

# p+ip 手征超导（Read-Green，D 类）
# H(k) = ξ_k σ_z + Δ(sin kx σ_x + sin ky σ_y)，ξ_k = -2t(cos kx+cos ky) - μ
t = 1.0
Delta = 0.5

def xi(kx, ky, mu):
    return -2*t*(np.cos(kx) + np.cos(ky)) - mu

def chern_number(mu, Ngrid=200):
    """算 p+ip 的陈数（下带 Berry 曲率积分）"""
    k = np.linspace(-np.pi, np.pi, Ngrid, endpoint=False)
    dk = 2*np.pi/Ngrid
    C = 0.0
    for kx in k:
        for ky in k:
            d = np.array([Delta*np.sin(kx), Delta*np.sin(ky), xi(kx, ky, mu)])
            norm = np.linalg.norm(d)
            dhat = d/norm
            # 有限差分算 Berry 曲率
            eps = 1e-4
            dxp = dhat_fd(kx+eps, ky, mu)
            dxm = dhat_fd(kx-eps, ky, mu)
            dyp = dhat_fd(kx, ky+eps, mu)
            dym = dhat_fd(kx, ky-eps, mu)
            F = np.dot(dhat, np.cross((dxp-dxm)/(2*eps), (dyp-dym)/(2*eps)))
            C += F * dk*dk
    return C / (2*np.pi)

def dhat_fd(kx, ky, mu):
    d = np.array([Delta*np.sin(kx), Delta*np.sin(ky), xi(kx, ky, mu)])
    return d/np.linalg.norm(d)

def carrier_density(mu, Ngrid=400):
    """正常态费米海体积（载流子密度 n）"""
    k = np.linspace(-np.pi, np.pi, Ngrid, endpoint=False)
    dk = 2*np.pi/Ngrid
    n = 0.0
    for kx in k:
        for ky in k:
            if xi(kx, ky, mu) < 0:
                n += dk*dk
    return n / (2*np.pi)**2   # 0~1 填充

mus = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
print("\n[扫 μ] 陈数 C 和 载流子密度 n：")
print("  μ      C(陈数)    n(填充)")
for mu in mus:
    C = chern_number(mu, Ngrid=60)
    n = carrier_density(mu, Ngrid=100)
    print(f"  {mu:+.0f}    {C:+.2f}      {n:.3f}")

print("\n>>> 陈数 C 是「拓扑阶跃」（|μ|<4t 时 C=1，否则 0）。")
print(">>> 载流子密度 n 是「光滑费米海」（从 0 平滑涨到 1）。")
print(">>> 两者【不成比例】：n 光滑连续，C 是 0/1 阶跃。")

# 关键结论
print("\n" + "="*72)
print("结论")
print("="*72)
print("标准 p+ip 里，「载流子密度 n」是费米海体积（光滑），")
print("「陈数 C」是拓扑不变量（阶跃）——两者独立、不成比例。")
print()
print("所以「物质=缺陷=拓扑」里的「n = 拓扑荷密度」在标准体模型里不成立：")
print("    n 不 ∝ C（n 光滑、C 阶跃）。")
print(">>> 「拓扑因子 f(C) ∝ C^2」不是从「载流子密度 n ∝ C」来的——")
print("    它需要一个【额外机制】（涡旋贡献 / 量子几何二次 / Majorana 杂化），")
print("    而标准 p+ip 体模型里没有这个机制。")

import numpy as np

print("="*72)
print("量子几何贡献 ∫g 对陈数 C 的标度：C、C² 还是模型依赖？")
print("="*72)

# Qi-Wu-Zhang 陈绝缘体：d(k) = (sin kx, sin ky, m + cos kx + cos ky)
# 量子度规 g_ij = (1/4) ∂i d̂ · ∂j d̂，超流权重（平带）D_s = (1/2)∫g
def dvec(kx, ky, m):
    return np.array([np.sin(kx), np.sin(ky), m + np.cos(kx) + np.cos(ky)])

def chern_and_metric(m, Ngrid=150):
    k = np.linspace(-np.pi, np.pi, Ngrid, endpoint=False)
    dk = 2*np.pi/Ngrid
    C = 0.0
    G = 0.0   # ∫ g_xx
    eps = 1e-4
    for kx in k:
        for ky in k:
            d = dvec(kx, ky, m)
            n = np.linalg.norm(d)
            if n < 1e-10: continue
            dh = d/n
            dx = (dvec(kx+eps,ky,m)-dvec(kx-eps,ky,m))/(2*eps)
            dy = (dvec(kx,ky+eps,m)-dvec(kx,ky-eps,m))/(2*eps)
            # 量子度规（去掉对 d̂ 方向的投影）
            dhdx = dx/n - dh*(np.dot(dh, dx/n))
            dhdy = dy/n - dh*(np.dot(dh, dy/n))
            gxx = np.dot(dhdx, dhdx)
            F = np.dot(dh, np.cross(dx/n, dy/n))
            C += F * dk*dk
            G += gxx * dk*dk
    return C/(4*np.pi), G/(4*(2*np.pi)**2)

print("\n  m      C(陈数)    ∫g(量子度规)   ∫g/C(比值)")
for m in [0.0, 0.5, 1.0, 1.5]:
    C, G = chern_and_metric(m)
    ratio = G/abs(C) if abs(C)>1e-6 else float('nan')
    print(f"  {m:.1f}    {C:+.2f}      {G:.4f}       {ratio:.3f}")

print("\n>>> 陈数 C 取 0/1（阶跃），量子度规 ∫g 随 m 平滑变化。")
print(">>> ∫g 和 C 不是「∝C」也不是「∝C²」——是「模型依赖」，且 ∫g ≥ |C|（下界）。")

# 关键：∫g 有没有 C² 的迹象？
print("\n[关键判定]")
print("  f(C)∝C² 需要「∫g 或 D_s 随 C 二次增长」。")
print("  但这里 C 只是 0/1 阶跃，∫g 是模型依赖的平滑量——给不出 C²。")
print("  要 C²，需要「大陈数 C≫1」且「∫g ∝ C²」的特定模型——非通用。")

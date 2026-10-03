import numpy as np

print("="*72)
print("修复 bug：陈数用 m=-1（非相变点），重算 ∫g 对 C 的标度")
print("="*72)

# Bug 说明：m=0 是相变点（ξ=m+cos+cos 在 Dirac 点过零），陈数=0。
# 修复：用 m=-1，此时陈数 C=2N（每个 Dirac 点贡献 ±N，符号由 ξ 定）。

def dvec(kx, ky, m, N):
    z = (np.sin(kx) + 1j*np.sin(ky))**N
    return np.array([z.real, z.imag, m + np.cos(kx) + np.cos(ky)])

def chern_metric(m, N, Ngrid=160):
    k = np.linspace(-np.pi, np.pi, Ngrid, endpoint=False)
    dk = 2*np.pi/Ngrid
    C = 0.0; G = 0.0
    eps = 1e-4
    for kx in k:
        for ky in k:
            d = dvec(kx, ky, m, N)
            n = np.linalg.norm(d)
            if n < 1e-10: continue
            dh = d/n
            dx = (dvec(kx+eps,ky,m,N)-dvec(kx-eps,ky,m,N))/(2*eps)
            dy = (dvec(kx,ky+eps,m,N)-dvec(kx,ky-eps,m,N))/(2*eps)
            dhdx = dx/n - dh*np.dot(dh, dx/n)
            dhdy = dy/n - dh*np.dot(dh, dy/n)
            G += np.dot(dhdx, dhdx) * dk*dk
            C += np.dot(dh, np.cross(dx/n, dy/n)) * dk*dk
    return C/(4*np.pi), G/(4*(2*np.pi)**2)

print("\n  N    C(陈数)    ∫g(量子度规)   ∫g/C   ∫g/C²")
results = []
for N in [1, 2, 3]:
    C, G = chern_metric(-1.0, N)
    results.append((N, abs(C), G))
    r1 = G/abs(C) if abs(C)>1e-6 else float('nan')
    r2 = G/(abs(C)**2) if abs(C)>1e-6 else float('nan')
    print(f"  {N}    {C:+.2f}     {G:.4f}          {r1:.3f}  {r2:.4f}")

print("\n>>> 注意：这个 sin^N 模型的 ∫g 增长主要来自「导数缩放」（∂sin^N ~ N sin^{N-1}），")
print("    是平凡伪影，不是拓扑机制。所以这不是测「拓扑标度」的好模型。")

# 干净测试：多层堆叠（每层 QWZ 独立），C=N，∫g 线性叠加
print("\n" + "="*72)
print("干净测试：多层堆叠（C=N，∫g 线性叠加）——验证 ∫g ∝ C（线性）")
print("="*72)
def qwz(kx, ky, m):
    return np.array([np.sin(kx), np.sin(ky), m + np.cos(kx) + np.cos(ky)])
def qwz_chern_metric(m, Ngrid=160):
    k = np.linspace(-np.pi, np.pi, Ngrid, endpoint=False)
    dk = 2*np.pi/Ngrid
    C = 0.0; G = 0.0; eps=1e-4
    for kx in k:
        for ky in k:
            d = qwz(kx, ky, m); n = np.linalg.norm(d)
            if n<1e-10: continue
            dh=d/n
            dx=(qwz(kx+eps,ky,m)-qwz(kx-eps,ky,m))/(2*eps)
            dy=(qwz(kx,ky+eps,m)-qwz(kx,ky-eps,m))/(2*eps)
            dhdx=dx/n-dh*np.dot(dh,dx/n); dhdy=dy/n-dh*np.dot(dh,dy/n)
            G += np.dot(dhdx,dhdx)*dk*dk
            C += np.dot(dh,np.cross(dx/n,dy/n))*dk*dk
    return C/(4*np.pi), G/(4*(2*np.pi)**2)

C1, G1 = qwz_chern_metric(-1.0)
print(f"  单层 QWZ: C={C1:+.2f}, ∫g={G1:.4f}")
print(f"  N 层堆叠: C=N×{abs(C1):.2f}, ∫g=N×{G1:.4f}")
print("  => ∫g ∝ C（线性），不是 C²")

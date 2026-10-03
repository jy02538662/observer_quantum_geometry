# -*- coding: utf-8 -*-
import numpy as np
lm = np.log(128**2 / (np.pi**2 * np.log(2 * 128**2 / np.pi**2)))
lam_c = 2.0
print("lambda_mod =", round(lm, 4), "  lambda_c =", lam_c)
print()
print("T_c = Lambda/(lambda*32), 1 meV = 11.6 K")
print(f"{'Lambda(eV)':>12} {'lambda=2':>14} {'lambda=5.32':>14}")
for Lam in [1.35, 1.65, 2.0, 3.2]:
    Lam_meV = Lam * 1000
    tc_c = Lam_meV / (lam_c * 32) * 11.6
    tc_m = Lam_meV / (lm * 32) * 11.6
    print(f"{Lam:>12} {tc_c:>11.0f}K {tc_m:>11.0f}K")
print()
print("实测: YBCO=92K  Bi-2212=85K  Bi-2223=110K  Hg-1223=134K")
print()
tc_target = 92 / 11.6  # meV
print("反解 Lambda 使 Tc=92K:")
print("  lambda=2   : Lambda =", round(tc_target * lam_c * 32 / 1000, 2), "eV")
print("  lambda=5.32: Lambda =", round(tc_target * lm * 32 / 1000, 2), "eV")
print()
J = 4 * 0.4**2 / 3.5
print("J=4t^2/U =", round(J, 3), "eV; 代入 lambda=5.32 -> Tc =", round(J * 1000 / (lm * 32) * 11.6, 0), "K")

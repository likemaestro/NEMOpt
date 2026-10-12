# Evidence for NDR-0065 (main strategy report): linear 2D constants, Hashin-Shtrikman bulk
# bound and p = 3 SIMP compatibility under plane strain vs plane stress. NumPy only.
# Run: python3 -I plane_stress_comparison.py
import numpy as np

E, nu = 100.0, 0.3
mu = E / (2 * (1 + nu))


def hs_upper(f, k1, m1, k2, m2):
    kH = k1 + (1 - f) / (1 / (k2 - k1) + f / (k1 + m1))
    z = m1 * (k1 + 2 * m1) / (2 * (k1 + m1))
    mH = m1 + (1 - f) / (1 / (m2 - m1) + f / (m1 + z))
    return kH, mH


cases = (("plane strain", E * nu / ((1 + nu) * (1 - 2 * nu)) + mu),   # 2D bulk = lambda + mu
         ("plane stress", E / (2 * (1 - nu))))                        # 2D bulk = E / (2(1 - nu))
rho = np.linspace(0.001, 0.999, 999)
for name, k in cases:
    nu2 = (k - mu) / (k + mu)                                         # in-plane Poisson's ratio
    ratio = max(r**3 * k / hs_upper(r, k, mu, 1e-12 * k, 1e-12 * mu)[0] for r in rho)
    print(f"{name:12s} k2D={k:.2f} nu_in-plane={nu2:.4f} "
          f"p_min=max(2/(1-nu),4/(1+nu))={max(2/(1-nu2), 4/(1+nu2)):.3f} "
          f"max SIMP(p=3)/HS bulk={ratio:.4f} HS+ k(f=0.4)={hs_upper(0.4, k, mu, 1e-12*k, 1e-12*mu)[0]:.3f}")

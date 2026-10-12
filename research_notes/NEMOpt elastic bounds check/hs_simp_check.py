# Evidence for NDR-0064 (main strategy report): 2D Hashin-Shtrikman upper bounds for the
# FEAD2026 constants and the SIMP/HS compatibility check. NumPy only.
# Run: python3 -I hs_simp_check.py
import numpy as np

E, nu = 100.0, 0.3                        # FEAD2026 solid, plane strain
lam = E * nu / ((1 + nu) * (1 - 2 * nu))
mu = E / (2 * (1 + nu))
k = lam + mu                              # 2D (plane-strain) bulk modulus
nu2d = lam / (lam + 2 * mu)               # in-plane Poisson's ratio
print(f"lambda={lam:.4f} mu={mu:.4f} k2D={k:.4f} nu2D={nu2d:.4f}")


def hs_upper(f, k1, m1, k2, m2):
    """2D Hashin-Shtrikman upper bounds (phase 1 stiffer, volume fraction f)."""
    kH = k1 + (1 - f) / (1 / (k2 - k1) + f / (k1 + m1))
    z = k1 * m1 / (k1 + 2 * m1)  # 2D Hashin-Shtrikman: zeta = kappa*mu/(kappa + 2mu) (corrected)
    mH = m1 + (1 - f) / (1 / (m2 - m1) + f / (m1 + z))
    return kH, mH


tiny, eps = 1e-12, 5e-4                   # void; FEAD2026 ersatz floor
for f in (0.4, 0.5, 0.6):
    kv, mv = hs_upper(f, k, mu, tiny * k, tiny * mu)
    ke, me = hs_upper(f, k, mu, eps * k, eps * mu)
    print(f"f={f:.1f}: Voigt k={f*k:.3f} | HS+ k void={kv:.4f} eps={ke:.4f} | HS+ mu void={mv:.4f} eps={me:.4f}")
print(f"kappa=28.4368 / HS+(0.6, void) = {28.4368 / hs_upper(0.6, k, mu, tiny*k, tiny*mu)[0]:.4f}")

print("Compatibility formula p >= max(2/(1-nu0), 4/(1+nu0)) (Bendsoe & Sigmund 1999, Eq. 11, p. 639; nu0 = in-plane Poisson ratio):")
print(f"  nu0 = 0.3 -> {max(2/(1-nu), 4/(1+nu)):.3f};  nu0 = nu2D -> {max(2/(1-nu2d), 4/(1+nu2d)):.3f}")
rho = np.linspace(0.001, 0.999, 999)
for p in (3.0, 3.5):
    rk = np.array([r**p * k / hs_upper(r, k, mu, tiny*k, tiny*mu)[0] for r in rho])
    rm = np.array([r**p * mu / hs_upper(r, k, mu, tiny*k, tiny*mu)[1] for r in rho])
    print(f"p={p}: max SIMP/HS bulk {rk.max():.4f} at rho={rho[rk.argmax()]:.3f}; shear {rm.max():.4f} at rho={rho[rm.argmax()]:.3f}")

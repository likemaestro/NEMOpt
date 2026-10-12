import numpy as np
E, nu = 100.0, 0.3
lam = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu)); k = lam + mu
def mu_hs_upper(f, k1, m1, m2, zeta):
    return m1 + (1-f)/(1/(m2-m1) + f/(m1+zeta))
z_script = mu*(k+2*mu)/(2*(k+mu))          # as in hs_simp_check.py
z_2d = k*mu/(k+2*mu)                       # 2D HS (Hashin 1965 transverse shear; Milton trace formula with d=2)
tiny = 1e-12
for f in (0.4, 0.5, 0.6):
    print(f"f={f}: mu_HS+ script-z={mu_hs_upper(f,k,mu,tiny*mu,z_script):.4f}  2D-zeta={mu_hs_upper(f,k,mu,tiny*mu,z_2d):.4f}")
# cross-check 2D zeta against trace form: 1/(m1+zeta) == (k1+2m1)/(2 m1 (k1+m1))
print("trace-form check:", 1/(mu+z_2d), (k+2*mu)/(2*mu*(k+mu)))
rho = np.linspace(0.001, 0.999, 999)
for p in (3.0, 3.5):
    for name, z in (("script", z_script), ("2D", z_2d)):
        r = np.array([x**p*mu/mu_hs_upper(x,k,mu,tiny*mu,z) for x in rho])
        print(f"p={p} {name}: max SIMP/HS shear {r.max():.4f} at rho={rho[r.argmax()]:.3f}")

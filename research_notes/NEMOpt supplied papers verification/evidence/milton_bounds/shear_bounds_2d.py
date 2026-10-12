# 2D shear-type bounds, FEAD2026 solid (plane strain E=100, nu=0.3), void or eps-material.
# (i) Milton Eq. 23.30 (optimal energy upper bound, any geometry) at pure shear eps0 eigenvalues (1,-1):
#     bounds (1/4) eps0:C*:eps0, e.g. C1212 (eps0 = e1e2+e2e1) or (C1111+C2222-2C1122)/4.
# (ii) 2D trace bound (Milton Eq. 23.33/23.44 pattern, d=2, M = Lambda_s, reference C1), derived:
#     f2*[1/(2(mu1-mu*_1)) + 1/(2(mu1-mu*_2))] <= 1/(mu1-mu2) - f1*(k1+2mu1)/(2mu1(k1+mu1))
# (iii) isotropic 2D HS shear upper bound.
E, nu = 100.0, 0.3
lam = E*nu/((1+nu)*(1-2*nu)); mu1 = E/(2*(1+nu)); k1 = lam + mu1
def e2330_pure_shear(f1, k1, mu1, k2, mu2):
    f2 = 1-f1
    return (mu1*mu2 + k1*(f1*mu1+f2*mu2))/(k1+f1*mu2+f2*mu1)   # third regime, s=0, d=2 -> /4
def hs_shear_2d(f1, k1, mu1, mu2):
    return mu1 + (1-f1)/(1/(mu2-mu1) + f1*(k1+2*mu1)/(2*mu1*(k1+mu1))) if mu2>0 else mu1 + (1-f1)/(-1/mu1 + f1*(k1+2*mu1)/(2*mu1*(k1+mu1)))
for eps in (0.0, 5e-4):
    for f1 in (0.4, 0.5, 0.6):
        k2, mu2 = eps*k1, eps*mu1
        rhs = (1/(mu1-mu2) - f1*(k1+2*mu1)/(2*mu1*(k1+mu1)))
        print(f"eps={eps:g} f1={f1}: single shear (23.30) <= {e2330_pure_shear(f1,k1,mu1,k2,mu2):.3f}; "
              f"isotropic HS shear <= {hs_shear_2d(f1,k1,mu1,mu2):.3f}; Voigt f1*mu1+f2*mu2 = {f1*mu1+(1-f1)*mu2:.3f}; "
              f"trace-bound RHS/f2 = {rhs/(1-f1):.6f}")

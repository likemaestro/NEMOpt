# Check Milton (2002) Eq. 23.30 (2D optimal upper energy bound, well-ordered phases)
# at hydrostatic average strain eps0 = I (eps1 = eps2 = 1), against the 2D HS upper bulk bound.
E, nu = 100.0, 0.3
lam = E*nu/((1+nu)*(1-2*nu)); mu1 = E/(2*(1+nu)); k1 = lam + mu1   # plane strain, 2D bulk modulus
def milton_2330_upper(f1, k1, mu1, k2, mu2, e1=1.0, e2=1.0):
    f2 = 1-f1; s, d = abs(e1+e2), abs(e1-e2)
    if (mu1 + f1*k2 + f2*k1)*d >= f1*(k1-k2)*s and (k1 + f1*mu2 + f2*mu1)*s >= f1*(mu1-mu2)*d:
        return s**2*(f1*k1+f2*k2) + d**2*(f1*mu1+f2*mu2) - f1*f2*((k1-k2)*s-(mu1-mu2)*d)**2/(f1*(mu2+k2)+f2*(mu1+k1))
    if (mu1 + f1*k2 + f2*k1)*d <= f1*(k1-k2)*s:
        return mu1*d**2 + (k1*k2 + mu1*(f1*k1+f2*k2))/(mu1+f1*k2+f2*k1)*s**2
    return k1*s**2 + (mu1*mu2 + k1*(f1*mu1+f2*mu2))/(k1+f1*mu2+f2*mu1)*d**2
def hs_upper_bulk_2d(f1, k1, mu1, k2):
    f2 = 1-f1
    return k1 + f2/(1/(k2-k1) + f1/(k1+mu1))
for eps in (0.0, 5e-4):
    for f1 in (0.4, 0.5, 0.6):
        k2, mu2 = eps*k1, eps*mu1
        W = milton_2330_upper(f1, k1, mu1, k2, mu2)   # = I:C*:I upper bound = 4*kappa
        print(f"eps={eps:g} f1={f1}: Milton 23.30 /4 = {W/4:.4f}   HS 2D upper = {hs_upper_bulk_2d(f1,k1,mu1,k2) if k2>0 else k1 + (1-f1)/(-1/k1 + f1/(k1+mu1)):.4f}")
print("k1, mu1 =", round(k1,4), round(mu1,4))

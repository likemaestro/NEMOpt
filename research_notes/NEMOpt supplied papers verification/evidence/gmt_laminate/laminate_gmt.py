# Reproduce GMT 1993 Sec. 6.5-6.6 macroscopic loss of ellipticity for a 2-phase laminate.
# Layers normal to x2, stretch lam along x1 (layers), lateral stress sigma22 = 0 in every layer.
# Energy (6.52): W = mu/2 (v1^2+v2^2-2-2 ln(v1 v2)) + kap/2 (v1 v2 - 1)^2, kap = nu E/((1+nu)(1-2nu)).
# Variant "bw": NEMOpt/Bonet-Wood volumetric term lam_L/2 (ln J)^2 instead of kap/2 (J-1)^2.
import math

def moduli(E, nu):
    mu = E / (2 * (1 + nu))
    kap = nu * E / ((1 + nu) * (1 - 2 * nu))
    return mu, kap

def sigma22(v1, v2, mu, kap, law):
    J = v1 * v2
    if law == "gmt":
        # (6.54): sigma22 = mu/(v1 v2) (v2^2 - 1) + kap (v1 v2 - 1)
        return mu / J * (v2 * v2 - 1) + kap * (J - 1)
    # Bonet-Wood: P = mu F - (mu - lamL ln J) F^-T ; sigma = P F^T / J
    return (mu * (v2 * v2 - 1) + kap * math.log(J)) / J

def solve_v2(v1, mu, kap, law):
    if law == "gmt":
        r = mu / kap
        return (v1 + math.sqrt(v1 * v1 + 4 * r * (v1 * v1 + r))) / (2 * (v1 * v1 + r))  # (6.55)
    lo, hi = 0.2, 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sigma22(v1, mid, mu, kap, law) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)

def L(v1, v2, mu, kap, law):
    J = v1 * v2
    if law == "gmt":  # (6.56)
        return dict(L1111=mu * (1 + 1 / v1**2) + kap * v2**2,
                    L1122=kap * (2 * v1 * v2 - 1),
                    L2222=mu * (1 + 1 / v2**2) + kap * v1**2,
                    L1212=mu, L2121=mu,
                    L1221=mu / J + kap * (1 - J))
    # Bonet-Wood nominal tangent A_iJkL = mu d_ik d_JL + lamL Finv_Ji Finv_Lk + (mu - lamL ln J) Finv_Li Finv_Jk
    f1, f2 = 1 / v1, 1 / v2
    c = mu - kap * math.log(J)
    return dict(L1111=mu + kap * f1 * f1 + c * f1 * f1,
                L1122=kap * f1 * f2,
                L2222=mu + kap * f2 * f2 + c * f2 * f2,
                L1212=mu, L2121=mu,
                L1221=c * f1 * f2)

def homogenize(layers):
    # (6.6); <f> = sum H_p f_p
    avg = lambda g: sum(H * g(Lp) for H, Lp in layers)
    iL22 = avg(lambda l: 1 / l["L2222"])
    iL12 = avg(lambda l: 1 / l["L1212"])
    a22 = avg(lambda l: l["L1122"] / l["L2222"])
    a12 = avg(lambda l: l["L1221"] / l["L1212"])
    return dict(
        M1111=avg(lambda l: l["L1111"] - l["L1122"] ** 2 / l["L2222"]) + a22 ** 2 / iL22,
        M1122=a22 / iL22,
        M2222=1 / iL22,
        M2121=avg(lambda l: l["L2121"] - l["L1221"] ** 2 / l["L1212"]) + a12 ** 2 / iL12,
        M1221=a12 / iL12,
        M1212=1 / iL12)

def alphaM(M, n=7200):
    # min over unit b of the min eigenvalue of Q_ik = M_iJkL b_J b_L (= (6.42) scan)
    best = 1e30
    for i in range(n + 1):
        p = 0.5 * math.pi * i / n
        b1, b2 = math.cos(p), math.sin(p)
        q11 = M["M1111"] * b1 * b1 + M["M1212"] * b2 * b2
        q22 = M["M2121"] * b1 * b1 + M["M2222"] * b2 * b2
        q12 = (M["M1122"] + M["M1221"]) * b1 * b2
        lmin = 0.5 * (q11 + q22) - math.sqrt(0.25 * (q11 - q22) ** 2 + q12 * q12)
        best = min(best, lmin)
    return best

def first_loss(Hb, law, E=1.0):
    soft = moduli(0.1 * E, 0.33)
    stiff = moduli(1.5 * E, 0.48)
    def am(lam):
        layers = []
        for H, (mu, kap) in ((1 - Hb, soft), (Hb, stiff)):
            v2 = solve_v2(lam, mu, kap, law)
            layers.append((H, L(lam, v2, mu, kap, law)))
        return alphaM(homogenize(layers), n=1800)
    lam, step = 1.0, 0.001
    prev = am(lam)
    while lam > 0.5:
        nl = lam - step
        cur = am(nl)
        if cur <= 0:
            lo, hi = nl, lam
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if am(mid) > 0:
                    hi = mid
                else:
                    lo = mid
            return 0.5 * (lo + hi)
        lam, prev = nl, cur
    return float("nan")

if __name__ == "__main__":
    # sanity: at lam = 1 the tangent equals isotropic linear elasticity
    mu, kap = moduli(1.0, 0.3)
    print("L at identity (gmt):", L(1, 1, mu, kap, "gmt"))
    for Hb in (0.1, 0.2, 0.3, 0.5, 0.7):
        print(f"H_b={Hb}: lambda_H^- (GMT law) = {first_loss(Hb, 'gmt'):.4f}; (Bonet-Wood volumetric) = {first_loss(Hb, 'bw'):.4f}")
    # incompressible-limit formula (1 - mu_harm/mu_arith)^(1/4) with GMT's shear moduli
    for Hb in (0.1, 0.5):
        ms, mf = moduli(0.1, 0.33)[0], moduli(1.5, 0.48)[0]
        mbar = Hb * mf + (1 - Hb) * ms
        mbrev = 1 / (Hb / mf + (1 - Hb) / ms)
        print(f"H_b={Hb}: incompressible formula with same mu ratio ({mf/ms:.3f}) gives {(1 - mbrev / mbar) ** 0.25:.4f}")

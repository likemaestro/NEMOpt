# Evidence for the Geymonat, Mueller & Triantafyllidis (1993) verification notes:
# macroscopic loss-of-ellipticity stretch of a two-phase laminate with NEMOpt-like phases.
# Run from this folder:  python3 -I laminate_nemopt.py laminate_gmt.py
# (laminate_gmt.py alone reproduces GMT's quoted lambda = 0.927 at stiff fraction 0.5).
import importlib.util, sys, math
spec = importlib.util.spec_from_file_location("lg", sys.argv[1]); lg = importlib.util.module_from_spec(spec); spec.loader.exec_module(lg)
def first_loss(Hb, ratio, nu, law):
    soft = lg.moduli(1.0, nu); stiff = lg.moduli(ratio, nu)
    def am(lam):
        layers = []
        for H, (mu, kap) in ((1 - Hb, soft), (Hb, stiff)):
            v2 = lg.solve_v2(lam, mu, kap, law)
            layers.append((H, lg.L(lam, v2, mu, kap, law)))
        return lg.alphaM(lg.homogenize(layers), n=1800)
    lam = 1.0
    while lam > 0.5:
        nl = lam - 0.001
        if am(nl) <= 0:
            lo, hi = nl, lam
            for _ in range(40):
                mid = 0.5*(lo+hi)
                if am(mid) > 0: hi = mid
                else: lo = mid
            return 0.5*(lo+hi)
        lam = nl
    return float('nan')
for nu in (0.3, 0.45, 0.49, 0.499):
    for Hb in (0.5, 0.2):
        print(f"nu={nu} ratio=10 vf={Hb}: BW law {first_loss(Hb,10,nu,'bw'):.4f}  GMT law {first_loss(Hb,10,nu,'gmt'):.4f}")
for Hb in (0.5, 0.2):
    mbar = Hb*10 + (1-Hb); mbrev = 1/(Hb/10 + (1-Hb))
    print(f"incompressible formula vf={Hb}: {(1-mbrev/mbar)**0.25:.4f}")

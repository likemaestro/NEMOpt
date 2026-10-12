import sys, numpy as np
sys.argv = ["x", "wt_mu8"]
src = open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "wt_bench.py")).read().split("np.set_printoptions")[0]
exec(src)
np.set_printoptions(precision=6, suppress=True, linewidth=140)
FM = np.array([[1.3,0],[0,0.8]]); groups = {1:[1,4,13,16], 2:[2,3,14,15], 3:[5,8,9,12], 4:[6,7,10,11]}
base = dens.copy(); eps = 1e-5
for gi, lst in groups.items():
    out = []
    for sgn in (+1, -1):
        dens[:] = base
        for e in lst: dens[e-1] += sgn*eps
        u, K, S = equilibrium(FM); AM = homog_A(FM, u, K)
        out.append((np.array([S[0,0],S[1,0],S[0,1],S[1,1]]), AM))
    dens[:] = base
    print("nu%d dS/dnu =" % gi, (out[0][0]-out[1][0])/(2*eps))
    print("dA/dnu =\n", (out[0][1]-out[1][1])/(2*eps))

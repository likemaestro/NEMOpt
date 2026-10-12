# Adjoint sensitivity (W&T Eq. 50-51) vs one-sided form (article Eq. 20) on W&T Sec. 4.1 case
import sys, numpy as np
sys.argv = ["x", "wt_mu8"]
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "wt_bench.py")).read().split("np.set_printoptions")[0])
def dtan(F, hh=1e-6):
    Z = np.zeros((2,2,2,2,2,2))
    for m in range(2):
        for n_ in range(2):
            Fp = F.copy(); Fp[m,n_] += hh; Fm = F.copy(); Fm[m,n_] -= hh
            Z[:,:,:,:,m,n_] = (tangent(Fp)-tangent(Fm))/(2*hh)
    return Z
FM = np.array([[1.3,0],[0,0.8]]); u, K, S = equilibrium(FM)
Es = [(0,0),(1,0),(0,1),(1,1)]
chis = {}
for (i,j) in Es:
    f = np.zeros(ndof)
    for e in range(n*n):
        d = edofs(e); ue = u[d].reshape(4,2)
        for xi in gp:
            for eta in gp:
                dNdX, dJ = grads(e, xi, eta); F = FM + ue.T @ dNdX; At = dens[e]**3*tangent(F)
                for a in range(4):
                    for k in range(2): f[d[2*a+k]] += (At[k,:,i,j] @ dNdX[a])*dJ
    chis[(i,j)] = solve_red(K, -f)
def G_at(e, dNdX, c):
    i,j = Es[c]; Eij = np.zeros((2,2)); Eij[i,j] = 1
    return Eij, Eij + chis[(i,j)][edofs(e)].reshape(4,2).T @ dNdX
groups = {1:[1,4,13,16], 2:[2,3,14,15], 3:[5,8,9,12], 4:[6,7,10,11]}
pairs = [(0,0),(3,3),(0,3),(1,1),(2,2),(1,2)]
np.set_printoptions(precision=6, suppress=True, linewidth=140)
for c1, c2 in pairs:
    f = np.zeros(ndof)
    for e in range(n*n):
        d = edofs(e); ue = u[d].reshape(4,2)
        for xi in gp:
            for eta in gp:
                dNdX, dJ = grads(e, xi, eta); F = FM + ue.T @ dNdX
                Z = dens[e]**3*dtan(F)
                _, G1 = G_at(e, dNdX, c1); _, G2 = G_at(e, dNdX, c2)
                Q = np.einsum('ijklmn,kl,mn->ij', Z, G1, G2)
                for a in range(4):
                    for k in range(2): f[d[2*a+k]] += (Q[k,:] @ dNdX[a])*dJ
    uh = solve_red(K, -f)
    sym = []; one = []
    for gi, lst in groups.items():
        s_sym = s_one = 0.0; vol = 4.0
        for el in lst:
            e = el-1; d = edofs(e); ue = u[d].reshape(4,2)
            for xi in gp:
                for eta in gp:
                    dNdX, dJ = grads(e, xi, eta); F = FM + ue.T @ dNdX
                    dA = 3*dens[e]**2*tangent(F); dS = 3*dens[e]**2*P_of_F(F)
                    E1, G1 = G_at(e, dNdX, c1); _, G2 = G_at(e, dNdX, c2)
                    gu = uh[d].reshape(4,2).T @ dNdX
                    s_sym += (np.einsum('ij,ijkl,kl->', G1, dA, G2) + np.sum(dS*gu))*dJ
                    s_one += (np.einsum('ij,ijkl,kl->', E1, dA, G2) + np.sum(dS*gu))*dJ
        sym.append(s_sym/vol); one.append(s_one/vol)
    lab = ["11","21","12","22"]
    print("dA_%s%s/dnu_i  symmetric (W&T Eq.51):" % (lab[c1], lab[c2]), np.array(sym), " one-sided (article Eq.20):", np.array(one))

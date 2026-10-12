# Reproduce Wallin & Tortorelli 2020 Sec. 4.1 (4x4 Q4 periodic cell, SIMP nu^3, no filter)
import numpy as np, sys
E, nu_m = 100.0, 0.3
lam = nu_m*E/((1+nu_m)*(1-2*nu_m)); mu = E/(2*(1+nu_m))
MODEL = sys.argv[1] if len(sys.argv) > 1 else "wt_mu8"
def P_of_F(F):
    J = F[0,0]*F[1,1]-F[0,1]*F[1,0]
    FinvT = np.array([[F[1,1], -F[1,0]], [-F[0,1], F[0,0]]])/J
    if MODEL.startswith("wt"):
        A = (lam+mu)/4; B = mu/8 if MODEL == "wt_mu8" else mu/2
        trC = np.sum(F*F)
        return A*(2*J*J-2)*FinvT + B*(2*F/J - trC/J*FinvT)
    else:  # NEMOpt / article Eq. 24
        return mu*F - (mu - lam*np.log(J))*FinvT
def tangent(F):
    h = 1e-30; Atan = np.zeros((2,2,2,2))
    for k in range(2):
        for l in range(2):
            Fc = F.astype(complex); Fc[k,l] += 1j*h
            Atan[:,:,k,l] = np.imag(P_of_F(Fc))/h
    return Atan
n = 4; L = 2.0; h = L/n
# nodes column-major from bottom-left as in W&T Fig. 3
xy = np.array([[-1+i*h, -1+j*h] for i in range(n+1) for j in range(n+1)])
nid = lambda i, j: i*(n+1)+j
els = []; dens = np.zeros(n*n)
for i in range(n):
    for j in range(n):
        els.append([nid(i,j), nid(i+1,j), nid(i+1,j+1), nid(i,j+1)])
val = {0.3:[1,4,13,16], 0.5:[2,3,14,15], 0.7:[5,8,9,12], 0.9:[6,7,10,11]}
for v, lst in val.items():
    for e in lst: dens[e-1] = v
# periodic map: node -> independent node
rep = {nid(i,j): nid(i % n, j % n) for i in range(n+1) for j in range(n+1)}
indep = sorted(set(rep.values())); idx = {nd:k for k,nd in enumerate(indep)}
fixed = idx[nid(2,2)]  # center node (0,0)
ndof = 2*len(indep)
gp = np.array([-1,1])/np.sqrt(3)
def shape(xi, eta):
    dN = 0.25*np.array([[-(1-eta), -(1-xi)], [(1-eta), -(1+xi)], [(1+eta), (1+xi)], [-(1+eta), (1-xi)]])
    return dN
def edofs(e):
    d = []
    for nd in els[e]:
        k = idx[rep[nd]]; d += [2*k, 2*k+1]
    return d
def grads(e, xi, eta):
    dN = shape(xi, eta); Jm = xy[els[e]].T @ dN; dNdX = dN @ np.linalg.inv(Jm)
    return dNdX, np.linalg.det(Jm)
def assemble(FM, u, need_K=True):
    R = np.zeros(ndof); K = np.zeros((ndof, ndof)); Sbar = np.zeros((2,2)); vol = 0
    for e in range(n*n):
        d = edofs(e); ue = u[d].reshape(4,2)
        for xi in gp:
            for eta in gp:
                dNdX, dJ = grads(e, xi, eta)
                F = FM + ue.T @ dNdX
                s = dens[e]**3
                P = s*P_of_F(F); Sbar += P*dJ; vol += dJ
                # B: dF_ij/du_(a,k) = delta_ik dNdX[a,j]
                for a in range(4):
                    for k in range(2):
                        R[d[2*a+k]] += (P[k,:] @ dNdX[a])*dJ
                if need_K:
                    At = s*tangent(F)
                    for a in range(4):
                        for k in range(2):
                            for b in range(4):
                                for m in range(2):
                                    K[d[2*a+k], d[2*b+m]] += np.einsum('j,l,jl->', dNdX[a], dNdX[b], At[k,:,m,:])*dJ
    return R, K, Sbar/vol
def solve_red(K, rhs):
    keep = [i for i in range(ndof) if i not in (2*fixed, 2*fixed+1)]
    x = np.zeros(ndof); x[keep] = np.linalg.solve(K[np.ix_(keep,keep)], rhs[keep]); return x
def equilibrium(FM, nsteps=10):
    u = np.zeros(ndof)
    for t in np.linspace(0,1,nsteps+1)[1:]:
        F_t = np.eye(2) + t*(FM-np.eye(2))
        for it in range(50):
            R, K, _ = assemble(F_t, u)
            if np.linalg.norm(R) < 1e-12: break
            u -= solve_red(K, R)
    R, K, S = assemble(FM, u)
    return u, K, S
def homog_A(FM, u, K):
    # A_M[ij,kl] = <(E_ij+grad chi_ij) : A : (E_kl + grad chi_kl)>
    chis = {}
    Es = [(0,0),(1,0),(0,1),(1,1)]  # W&T order 11,21,12,22
    for (i,j) in Es:
        f = np.zeros(ndof)
        for e in range(n*n):
            d = edofs(e); ue = u[d].reshape(4,2)
            for xi in gp:
                for eta in gp:
                    dNdX, dJ = grads(e, xi, eta); F = FM + ue.T @ dNdX
                    At = dens[e]**3*tangent(F)
                    for a in range(4):
                        for k in range(2):
                            f[d[2*a+k]] += (At[k,:,i,j] @ dNdX[a])*dJ
        chis[(i,j)] = solve_red(K, -f)
    AM = np.zeros((4,4)); vol = 0
    for e in range(n*n):
        d = edofs(e); ue = u[d].reshape(4,2)
        for xi in gp:
            for eta in gp:
                dNdX, dJ = grads(e, xi, eta); F = FM + ue.T @ dNdX
                At = dens[e]**3*tangent(F); vol += dJ
                G = []
                for (i,j) in Es:
                    Eij = np.zeros((2,2)); Eij[i,j] = 1
                    G.append(Eij + chis[(i,j)][d].reshape(4,2).T @ dNdX)
                for p_ in range(4):
                    for q in range(4):
                        AM[p_,q] += np.einsum('ij,ijkl,kl->', G[p_], At, G[q])*dJ
    return AM/vol
np.set_printoptions(precision=5, suppress=True, linewidth=120)
FM = np.array([[1.3,0],[0,0.8]])
u, K, S = equilibrium(FM)
print(MODEL, "S_M [S11 S21 S12 S22] =", np.array([S[0,0],S[1,0],S[0,1],S[1,1]]))
print("A_M (rows 11,21,12,22) =\n", homog_A(FM, u, K))
# uniaxial: find F22 with S22 = 0 (secant + Newton on scalar)
f22 = 0.85
for it in range(30):
    FM = np.array([[1.3,0],[0,f22]]); u, K, S = equilibrium(FM)
    AM = homog_A(FM, u, K)
    if abs(S[1,1]) < 1e-11: break
    f22 -= S[1,1]/AM[3,3]
print("uniaxial F22 =", f22, " S11 =", S[0,0], " S22 =", S[1,1])
print("A_M uniaxial =\n", AM)
CM = np.linalg.inv(AM); print("C_M =\n", CM)
print("E1 = 1/C1111 =", 1/CM[0,0], " nu12 = -E1*C2211 =", -CM[3,0]/CM[0,0])

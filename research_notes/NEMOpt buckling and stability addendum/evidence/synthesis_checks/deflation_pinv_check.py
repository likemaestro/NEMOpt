# Translation deflation without densifying: with G = zero-padded inverse of the center-pinned periodic tangent and
# P = I - T(T'T)^-1 T' (T = the two rigid translations), P*G*P equals the Moore-Penrose inverse of the unpinned
# singular tangent, so shift-invert Lanczos on x -> P*(G*(P*x)) returns the true nonzero spectrum without the two
# pinned-translation pseudo-modes. Heterogeneous n=10 cell. Run: python3 -I deflation_pinv_check.py
import numpy as np
exec(open('parity_inertia_check.py').read().split('# C2 action')[0])   # reuse assembly of K (unpinned, heterogeneous)
T=np.zeros((N,2)); T[0::2,0]=1; T[1::2,1]=1
P=np.eye(N)-T@np.linalg.solve(T.T@T,T.T)
c=node(0,0); pin=[2*c,2*c+1]; keep=[d for d in range(N) if d not in pin]
G=np.zeros((N,N)); G[np.ix_(keep,keep)]=np.linalg.inv(K[np.ix_(keep,keep)])
Kp=np.linalg.pinv(K,rcond=1e-12)
print("max|PGP - pinv(K)| / max|pinv(K)| =",np.abs(P@G@P-Kp).max()/np.abs(Kp).max())
ev=np.linalg.eigvalsh(K); evp=np.linalg.eigvalsh(K[np.ix_(keep,keep)])
mu=np.sort(1/np.linalg.eigvalsh(P@G@P)[-6:])
print("3 smallest nonzero eig of unpinned K :",ev[2:5])
print("3 smallest eig of pinned K           :",evp[:3])
print("3 smallest from 1/eig(PGP)           :",mu[:3])

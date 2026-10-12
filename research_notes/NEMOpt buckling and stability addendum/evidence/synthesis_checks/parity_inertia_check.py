# C2 parity split at k=0: inertia(full, center pinned) = inertia(odd sector) + inertia(even sector, center pinned)
# on a heterogeneous C2-symmetric linear-elastic periodic cell with shifts (surrogate for an indefinite tangent). Run: python3 -I parity_inertia_check.py
import numpy as np
class sl:
    @staticmethod
    def null_space(A,rcond=1e-10):
        u,sv,vh=np.linalg.svd(A)
        tol=rcond*max(sv.max() if sv.size else 0,1)
        r=int((sv>tol).sum())
        return vh[r:].conj().T

def q4ke(E,nu,h):
    lam=E*nu/((1+nu)*(1-2*nu)); mu=E/(2*(1+nu))
    C=np.array([[lam+2*mu,lam,0],[lam,lam+2*mu,0],[0,0,mu]])
    Ke=np.zeros((8,8)); g=1/np.sqrt(3)
    for xi,eta in ((-g,-g),(g,-g),(g,g),(-g,g)):
        dNdxi=0.25*np.array([-(1-eta),(1-eta),(1+eta),-(1+eta)])
        dNdeta=0.25*np.array([-(1-xi),-(1+xi),(1+xi),(1-xi)])
        J=h/2; dNdx=dNdxi/J; dNdy=dNdeta/J
        B=np.zeros((3,8))
        for a in range(4):
            B[0,2*a]=dNdx[a]; B[1,2*a+1]=dNdy[a]; B[2,2*a]=dNdy[a]; B[2,2*a+1]=dNdx[a]
        Ke+=B.T@C@B*J**2
    return Ke
n=10; h=1/n; Ke=q4ke(100.,0.3,h); N=2*n*n
rng=np.random.default_rng(5); r=0.02+0.98*rng.random((n,n))
rho=np.empty((n,n))
for i in range(n):
    for j in range(n):
        rho[i,j]=0.5*(r[i,j]+r[(-i-1)%n,(-j-1)%n])   # C2-symmetric (not mirror-symmetric) density
node=lambda i,j:(i%n)+n*(j%n)
K=np.zeros((N,N))
for j in range(n):
    for i in range(n):
        dofs=[]
        for a,b in ((i,j),(i+1,j),(i+1,j+1),(i,j+1)):
            nd=node(a,b); dofs+=[2*nd,2*nd+1]
        K[np.ix_(dofs,dofs)]+=rho[i,j]**3*Ke
# C2 action on nodal DOFs: (Rw)(y) = -w(-y)
R=np.zeros((N,N))
for j in range(n):
    for i in range(n):
        a=node(i,j); b=node(-i,-j)
        R[2*a,2*b]=-1; R[2*a+1,2*b+1]=-1
print("||KR-RK|| =",np.abs(K@R-R@K).max())
c=node(0,0); pin=[2*c,2*c+1]; keep=[d for d in range(N) if d not in pin]
ev=np.linalg.eigvalsh(K); 
odd=sl.null_space(R-np.eye(N)); even=sl.null_space(R+np.eye(N))   # odd: w(-y)=-w(y) (R=+1); even: R=-1
# pinned even sector: even fields with w(center)=0
E2=even[[d for d in range(N) if d not in pin],:]
Ep=sl.null_space(even[pin,:]); evenp=even@Ep
mism=0; tests=0
for sigma in np.quantile(ev[3:],[0.002,0.01,0.05,0.1,0.3,0.6]):
    Ks=K-sigma*np.eye(N)
    full=np.sum(np.linalg.eigvalsh(Ks[np.ix_(keep,keep)])<0)
    no=np.sum(np.linalg.eigvalsh(odd.T@Ks@odd)<0); ne=np.sum(np.linalg.eigvalsh(evenp.T@Ks@evenp)<0)
    tests+=1; mism+= (full!=no+ne)
    print(f"sigma={sigma:.4g}: neg(full pinned)={full}  neg(odd)={no}  neg(even, center pinned)={ne}  sum={no+ne}")
print("dims: odd",odd.shape[1],"even",even.shape[1],"even pinned",evenp.shape[1],"full pinned",len(keep))
print("mismatches:",mism,"of",tests)
# where does the lowest nonzero mode live?
w,V=np.linalg.eigh(K[np.ix_(keep,keep)])

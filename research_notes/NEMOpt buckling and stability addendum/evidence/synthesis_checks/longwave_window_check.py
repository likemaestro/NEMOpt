# Long-wave Bloch limit on a homogeneous n=12 periodic cell: lambda*n^2/theta^2 -> eigenvalues of Q(e1) (mu, lam+2mu).
# Also reproduces the pinned-vs-unpinned spectrum at theta=(pi,0). Run: python3 -I longwave_window_check.py
import numpy as np
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
    return Ke,lam,mu
def bloch(n,k,pin=False):
    h=1/n; Ke,lam,mu=q4ke(100.0,0.3,h); N=2*n*n
    K=np.zeros((N,N),dtype=complex)
    node=lambda i,j:(i%n)+n*(j%n)
    for j in range(n):
        for i in range(n):
            dofs=[];ph=[]
            for a,b in ((i,j),(i+1,j),(i+1,j+1),(i,j+1)):
                nd=node(a,b); s=np.exp(1j*(k[0]*(a>=n)+k[1]*(b>=n)))
                dofs+=[2*nd,2*nd+1]; ph+=[s,s]
            ph=np.array(ph)
            K[np.ix_(dofs,dofs)]+=np.conj(ph)[:,None]*Ke*ph[None,:]
    K=(K+K.conj().T)/2
    if pin:
        c=node(n//2,n//2); keep=[d for d in range(N) if d not in (2*c,2*c+1)]; K=K[np.ix_(keep,keep)]
    return K,lam,mu
n=12
for th in (1e-6,1e-4,1e-3,1e-2,1e-1):
    K,lam,mu=bloch(n,(th,0.0))
    ev=np.linalg.eigvalsh(K)[:2]
    print(f"theta={th:g}: two smallest eig = {ev[0]:.4e}, {ev[1]:.4e}; scaled by n^2/theta^2 = {ev[0]*n*n/th**2:.4f}, {ev[1]*n*n/th**2:.4f}  (continuum mu={mu:.4f}, lam+2mu={lam+2*mu:.4f})")
K,lam,mu=bloch(n,(np.pi,0.0)); print("unpinned (pi,0):",np.linalg.eigvalsh(K)[:3])
K,lam,mu=bloch(n,(np.pi,0.0),pin=True); print("pinned (pi,0):",np.linalg.eigvalsh(K)[:3])

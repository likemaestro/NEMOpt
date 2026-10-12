# Re-creates the exp1.jl n=12 pinned homogeneous periodic Q4 matrix (NumPy only) and checks
# whether the exp1 shifts sigma=(ev[m]+ev[m+1])/2, m=1,3,7, land on repeated eigenvalues. Run: python3 -I exp1_degeneracy_check.py
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
    return Ke
n=12; h=1/n; Ke=q4ke(100.0,0.3,h); N=2*n*n
K=np.zeros((N,N))
node=lambda i,j: (i%n)+n*(j%n)
for j in range(n):
    for i in range(n):
        dofs=[]
        for a,b in ((i,j),(i+1,j),(i+1,j+1),(i,j+1)):
            nd=node(a,b); dofs+= [2*nd,2*nd+1]
        K[np.ix_(dofs,dofs)]+=Ke
c=node(n//2,n//2); keep=[d for d in range(N) if d not in (2*c,2*c+1)]
K=K[np.ix_(keep,keep)]; K=(K+K.T)/2
ev=np.linalg.eigvalsh(K)
print("size",K.shape[0],"lambda_min",ev[0])
print("ev[1..10] (1-based):",["%.12g"%v for v in ev[:10]])
for m in (1,3,7):
    s=(ev[m-1]+ev[m])/2
    below=np.sum(ev< s - 1e-10*abs(s)); near=np.sum(np.abs(ev-s)<=1e-10*abs(s))
    print(f"m={m}: sigma={s:.12g}; ev[m]={ev[m-1]:.12g} ev[m+1]={ev[m]:.12g}; strictly below={below}, within 1e-10 rel of sigma={near}")

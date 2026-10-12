# FEAD2026 energy: closed-form min eig of the acoustic tensor vs angle scan; blended (H_Steep) law ellipticity. Run: python3 -I acoustic_closed_form_check.py
import numpy as np
E,nu=100.0,0.3; lam=E*nu/((1+nu)*(1-2*nu)); mu=E/(2*(1+nu))
d=np.eye(2)
def A_nh(F):
    Fi=np.linalg.inv(F).T; J=np.linalg.det(F); A=np.zeros((2,2,2,2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    A[i,j,k,l]=mu*d[i,k]*d[j,l]+lam*Fi[i,j]*Fi[k,l]+(mu-lam*np.log(J))*Fi[i,l]*Fi[k,j]
    return A
def C_lin():
    A=np.zeros((2,2,2,2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    A[i,j,k,l]=lam*d[i,j]*d[k,l]+mu*(d[i,k]*d[j,l]+d[i,l]*d[j,k])
    return A
def minQ(A,nt=7200):
    best=np.inf
    for t in np.linspace(0,np.pi,nt,endpoint=False):
        n=np.array([np.cos(t),np.sin(t)]); Q=np.einsum('ijkl,j,l->ik',A,n,n)
        best=min(best,np.linalg.eigvalsh((Q+Q.T)/2)[0])
    return best
def closed(F):
    J=np.linalg.det(F); s=np.linalg.svd(F,compute_uv=False); c=lam+mu-lam*np.log(J)
    return min(mu,mu+c/s.min()**2,mu+c/s.max()**2)
rng=np.random.default_rng(1); worst=0
for _ in range(200):
    F=np.eye(2)+0.4*rng.standard_normal((2,2))
    if np.linalg.det(F)<=0.05: continue
    a=minQ(A_nh(F),1440); b=closed(F); worst=max(worst,abs(a-b)/abs(b))
print("random F (200 draws): max rel diff scan(0.125deg) vs closed form =",worst)
# published load states
for lbl,F in (("1.075I",1.075*np.eye(2)),("case2 [1.025 .025;.025 1.025]",np.array([[1.025,.025],[.025,1.025]])),("bulk2 [1 .05;.05 1]",np.array([[1,.05],[.05,1.0]])),("1.0001I",1.0001*np.eye(2)),("uniax comp 0.9",np.diag([0.9,1.0]))):
    print(lbl,"closed-form min eig Q =",round(closed(F),6),"(mu =",round(mu,6),")")
# blended tangent with H_Steep weight w in [0,1] and g scale: positive combination -> elliptic
Cl=C_lin(); mins=[]
for _ in range(200):
    F=np.eye(2)+0.3*rng.standard_normal((2,2))
    if np.linalg.det(F)<=0.2: continue
    w=rng.uniform(); g=10**rng.uniform(-4,0)
    mins.append(minQ(g*(w*A_nh(F)+(1-w)*Cl),720)/g)
print("blended law, 200 random (F, w, g): min over draws of minQ/g =",min(mins))
print("uniaxial-tension LoE stretch exp(1+2mu/lam) =",np.exp(1+2*mu/lam))

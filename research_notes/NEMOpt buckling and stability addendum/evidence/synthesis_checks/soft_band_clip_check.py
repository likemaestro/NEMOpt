# Soft band of the FEAD2026 interpolation: filtered densities where the blend H_Steep(rho~^p) is nonlinear
# while the stiffness is clipped at g = eps; size of the clip kink and of smooth-floor void sensitivities.
# Run: python3 -I soft_band_clip_check.py
import numpy as np
def H(r,b,w=0.5): return (np.tanh(b*w)+np.tanh(b*(r-w)))/(np.tanh(b*w)+np.tanh(b*(1-w)))
def Hd(r,b,w=0.5): return b/np.cosh(b*(r-w))**2/(np.tanh(b*w)+np.tanh(b*(1-w)))
p=3; b2=200; w2=0.05; thr=w2**(1/p)
for eps in (5e-4,1e-5):
    rs=np.linspace(0,1,1000001); cb=rs[np.argmax(H(rs,32)**p>eps)]
    print(f"eps={eps}, p=3, beta=32: clip boundary rho~={cb:.4f} (H={H(cb,32):.4f}); blend threshold rho~={thr:.4f}")
    for r in (0.368,0.38,0.39,0.40,0.43,0.46):
        print(f"   rho~={r}: H_Steep(rho~^p)={H(r**p,b2,w2):.4f}, H^p={H(r,32)**p:.3e}, g={max(H(r,32)**p,eps):.3e}")
eps=5e-4
for beta in np.arange(8,12.01,0.25):
    if H(thr,beta)**3<=eps: print("band first appears at beta ~",beta); break
for beta in (16,24,32):
    rs=np.linspace(0,1,400001); cb=rs[np.argmax(H(rs,beta)**3>eps)]
    print(f"beta={beta}: band rho~ in [{thr:.3f}, {cb:.3f}]")
for (pp,beta) in ((1,1),(1.5,1),(3,1),(3,8),(3,32)):
    rs=np.linspace(0,1,200001); Hp=H(rs,beta)**pp; act=rs[Hp<=eps]
    gs=(1-eps)*pp*H(rs,beta)**(pp-1)*Hd(rs,beta)
    rb=act.max(); jump=pp*H(rb,beta)**(pp-1)*Hd(rb,beta)
    print(f"p={pp}, beta={beta}: clip active rho~<={rb:.4f}; g' jump at boundary {jump:.3e} ({jump/gs.max():.2e} of max g'); smooth-floor g' at rho~=0.3: {gs[np.argmin(abs(rs-0.3))]:.2e}")

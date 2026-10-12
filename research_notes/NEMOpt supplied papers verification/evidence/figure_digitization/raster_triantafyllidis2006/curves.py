import sys, json
from PIL import Image
import numpy as np
# args: image xL xR yT yB e1L e1R e2T e2B [excl boxes as json list of [e1min,e1max,e2min,e2max]]
img=sys.argv[1]
xL,xR,yT,yB=map(float,sys.argv[2:6])
e1L,e1R,e2T,e2B=map(float,sys.argv[6:10])
excl=json.loads(sys.argv[10]) if len(sys.argv)>10 else []
im=np.asarray(Image.open(img).convert("RGB")).astype(int)
ys,xs=np.mgrid[int(yT)+4:int(yB)-3, int(xL)+4:int(xR)-3]
sub=im[int(yT)+4:int(yB)-3, int(xL)+4:int(xR)-3]
r,g,b=sub[...,0],sub[...,1],sub[...,2]
cls={
 "red":  (r>170)&(g<90)&(b<70)&(r-g>110),
 "green":(g>90)&(r<110)&(g-r>50)&(g-b>30),
 "blue": (b>100)&(r<120)&(g<120)&(b-r>50)&(b-g>50),
 "magenta":(r>170)&(b>90)&(g<100),
 "black":(r+g+b<150),
}
e1=e1L+(xs-xL)/(xR-xL)*(e1R-e1L)
e2=e2T+(ys-yT)/(yB-yT)*(e2B-e2T)
def ok(m):
    for bx in excl:
        m=m&~((e1>=bx[0])&(e1<=bx[1])&(e2>=bx[2])&(e2<=bx[3]))
    return m
tol=float(sys.argv[11]) if len(sys.argv)>11 else 0.0025
for name,m in cls.items():
    m=ok(m)
    n=int(m.sum())
    if n<50: continue
    E1=e1[m]; E2=e2[m]
    out={"n":n}
    # near diagonal e1=e2
    d=np.abs(E1-E2)<tol
    if d.any(): out["diag(e1=e2)"]=sorted(set(np.round(E1[d],3)))
    t=np.abs(E2-e2T)<tol*1.5
    if t.any(): out["at e2=%g"%e2T]=sorted(set(np.round(E1[t],3)))
    rr=np.abs(E1-e1R)<tol*1.5
    if rr.any(): out["at e1=%g"%e1R]=sorted(set(np.round(E2[rr],3)))
    print(name, out)

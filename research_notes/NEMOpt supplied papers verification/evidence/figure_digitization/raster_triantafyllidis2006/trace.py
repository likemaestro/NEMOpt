import sys
from PIL import Image
import numpy as np
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
x0,x1,xv0,xv1 = [float(v) for v in sys.argv[2:6]]
y0,y1,yv0,yv1 = [float(v) for v in sys.argv[6:10]]
phis = [float(v) for v in sys.argv[10].split(',')]
lam_a, lam_b = float(sys.argv[11]), float(sys.argv[12])
mask_boxes = [tuple(int(t) for t in b.split(':')) for b in sys.argv[13].split(',')] if len(sys.argv)>13 and sys.argv[13] else []
def px(e1): return x0 + (e1-xv0)*(x1-x0)/(xv1-xv0)
def py(e2): return y0 + (e2-yv0)*(y1-y0)/(yv1-yv0)
H,W,_=im.shape
for phi in phis:
    c,s=np.cos(np.radians(phi)),np.sin(np.radians(phi))
    lam=np.linspace(lam_a,lam_b,6000)
    X=px(lam*c); Y=py(lam*s)
    ok=(X>min(x0,x1)+3)&(X<max(x0,x1)-3)&(Y>min(y0,y1)+3)&(Y<max(y0,y1)-3)
    Xi=np.round(X).astype(int); Yi=np.round(Y).astype(int)
    for (a,b,c2,d) in mask_boxes:
        ok &= ~((Xi>=a)&(Xi<=b)&(Yi>=c2)&(Yi<=d))
    Xi=np.clip(Xi,0,W-1); Yi=np.clip(Yi,0,H-1)
    col=im[Yi,Xi]
    dark=(col.sum(axis=1)<300)&ok
    red=(col[:,0]>150)&(col[:,1]<110)&(col[:,2]<110)&ok
    green=(col[:,1]>140)&(col[:,0]<150)&(col[:,2]<150)&ok
    blue=(col[:,2]>150)&(col[:,0]<110)&(col[:,1]<150)&ok
    def runs(m):
        out=[];i=0;n=len(m)
        while i<n:
            if m[i]:
                j=i
                while j<n and m[j]: j+=1
                if lam[j-1]-lam[i] > 0.0: out.append(round((lam[i]+lam[j-1])/2,4))
                else: out.append(round(lam[i],4))
                i=j
            else: i+=1
        return out
    res=f'phi={phi:6.1f}  black={runs(dark)}  red={runs(red)}'
    g=runs(green); b=runs(blue)
    if g: res+=f'  green={g}'
    if b: res+=f'  blue={b}'
    print(res)

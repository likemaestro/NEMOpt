import sys
from PIL import Image
import numpy as np
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
x0,x1,xv0,xv1 = [float(v) for v in sys.argv[2:6]]
y0,y1,yv0,yv1 = [float(v) for v in sys.argv[6:10]]
phi = float(sys.argv[10])  # ray angle in degrees in (e1,e2) plane
lam_a, lam_b = float(sys.argv[11]), float(sys.argv[12])
offs = [float(v) for v in sys.argv[13].split(',')] if len(sys.argv)>13 else [0.0]
c,s = np.cos(np.radians(phi)), np.sin(np.radians(phi))
def px(e1): return x0 + (e1-xv0)*(x1-x0)/(xv1-xv0)
def py(e2): return y0 + (e2-yv0)*(y1-y0)/(yv1-yv0)
for off in offs:
    # offset perpendicular to the ray, in data units
    lam = np.linspace(lam_a, lam_b, 4000)
    e1 = lam*c - off*s; e2 = lam*s + off*c
    X = np.round(px(e1)).astype(int); Y = np.round(py(e2)).astype(int)
    col = im[Y,X]
    dark = col.sum(axis=1) < 300
    red = (col[:,0]>150)&(col[:,1]<110)&(col[:,2]<110)
    def runs(mask):
        out=[]; i=0; n=len(mask)
        while i<n:
            if mask[i]:
                j=i
                while j<n and mask[j]: j+=1
                out.append((round(lam[i],4), round(lam[j-1],4), round((lam[i]+lam[j-1])/2,4)))
                i=j
            else: i+=1
        return out
    print(f'phi={phi} off={off}: dark runs (lam_start, lam_end, mid):', runs(dark))
    print(f'   red runs:', runs(red))

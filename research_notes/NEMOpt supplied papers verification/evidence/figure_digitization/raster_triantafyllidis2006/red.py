import sys
from PIL import Image
import numpy as np
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
# frame (b): x0 px->xv0, x1 px->xv1 ; y0 px->yv0, y1 px->yv1
x0,x1,xv0,xv1 = [float(v) for v in sys.argv[2:6]]
y0,y1,yv0,yv1 = [float(v) for v in sys.argv[6:10]]
split = float(sys.argv[10])
legend = [int(v) for v in sys.argv[11:15]]  # xmin,xmax,ymin,ymax to exclude
R,G,B = im[:,:,0],im[:,:,1],im[:,:,2]
red = (R>150)&(G<110)&(B<110)
ys,xs = np.where(red)
m = (xs>min(x0,x1)+2)&(xs<max(x0,x1)-2)&(ys>min(y0,y1)+2)&(ys<max(y0,y1)-2)
m &= ~((xs>=legend[0])&(xs<=legend[1])&(ys>=legend[2])&(ys<=legend[3]))
xs,ys = xs[m],ys[m]
e1 = xv0 + (xs-x0)*(xv1-xv0)/(x1-x0)
e2 = yv0 + (ys-y0)*(yv1-yv0)/(y1-y0)
s = e1+e2
for name,sel in [('outer(s<split)', s<split), ('inner(s>=split)', s>=split)]:
    a,b = e1[sel],e2[sel]
    print(name, 'n=',sel.sum())
    d = a-b
    near = np.abs(d)<0.002
    print('  diagonal crossing e1=e2 ~', np.round(np.mean(a[near]),4), 'spread', np.round(np.std(a[near]),4), 'n',near.sum())
    # nearest to e2 = top (yv0) and e1 = right (xv1)
    top = b > yv0 - 0.004
    if top.sum(): print('  at e2~0: e1 =', np.round(np.mean(a[top]),4), 'range', np.round(a[top].min(),4), np.round(a[top].max(),4))
    right = a > xv1 - 0.004
    if right.sum(): print('  at e1~0: e2 =', np.round(np.mean(b[right]),4), 'range', np.round(b[right].min(),4), np.round(b[right].max(),4))
    # fit e1+e2 stats
    print('  e1+e2 mean', np.round(np.mean(a+b),4), 'min', np.round((a+b).min(),4), 'max', np.round((a+b).max(),4))

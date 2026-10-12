import sys
from PIL import Image
import numpy as np
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
H,W,_ = im.shape
red = (im[:,:,0]>180)&(im[:,:,1]<120)&(im[:,:,2]<120)
green = (im[:,:,1]>150)&(im[:,:,0]<120)&(im[:,:,2]<120)
print(H,W)
# bounding box of red
ys,xs = np.where(red)
print('red bbox', xs.min(), xs.max(), ys.min(), ys.max())
def runs(v):
    idx = np.where(v)[0]
    if len(idx)==0: return []
    groups=[[idx[0]]]
    for i in idx[1:]:
        if i-groups[-1][-1]<=2: groups[-1].append(i)
        else: groups.append([i])
    return [int(np.mean(g)) for g in groups]
for y in [ys.max()-3, ys.max()-6, ys.min()+4, ys.min()+7]:
    r = runs(red[y, xs.min()+3:xs.max()-2])
    print('row',y,len(r))
for x in [xs.min()+4, xs.max()-4, xs.min()+7]:
    r = runs(red[ys.min()+3:ys.max()-2, x])
    print('col',x,len(r))
# green across center row/col
cy=(ys.min()+ys.max())//2; cx=(xs.min()+xs.max())//2
gy,gx=np.where(green)
print('green bbox', gx.min(), gx.max(), gy.min(), gy.max())
for y in [cy+5, cy-5]:
    print('green row', y, len(runs(green[y,:])))
for x in [cx+5]:
    print('green col', x, len(runs(green[:,x])))

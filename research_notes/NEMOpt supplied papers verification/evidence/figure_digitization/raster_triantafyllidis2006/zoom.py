import sys
from PIL import Image, ImageDraw
import numpy as np
src,out = sys.argv[1], sys.argv[2]
x0,x1,xv0,xv1 = [float(v) for v in sys.argv[3:7]]
y0,y1,yv0,yv1 = [float(v) for v in sys.argv[7:11]]
e1a,e1b,e2a,e2b,step,scale = [float(v) for v in sys.argv[11:17]]
im = Image.open(src).convert('RGB')
def px(e1): return x0 + (e1-xv0)*(x1-x0)/(xv1-xv0)
def py(e2): return y0 + (e2-yv0)*(y1-y0)/(yv1-yv0)
L,R = sorted([px(e1a),px(e1b)]); T,Bt = sorted([py(e2a),py(e2b)])
crop = im.crop((int(L),int(T),int(R)+1,int(Bt)+1))
crop = crop.resize((int(crop.width*scale), int(crop.height*scale)), Image.NEAREST)
d = ImageDraw.Draw(crop)
v = np.arange(np.ceil(min(e1a,e1b)/step)*step, max(e1a,e1b)+1e-12, step)
for e in v:
    X = (px(e)-int(L))*scale
    d.line([(X,0),(X,crop.height)], fill=(0,160,255), width=1)
    d.text((X+2,2), f"{e:.3f}", fill=(0,100,255))
v = np.arange(np.ceil(min(e2a,e2b)/step)*step, max(e2a,e2b)+1e-12, step)
for e in v:
    Y = (py(e)-int(T))*scale
    d.line([(0,Y),(crop.width,Y)], fill=(0,200,100), width=1)
    d.text((2,Y+2), f"{e:.3f}", fill=(0,150,80))
# diagonal
d.line([((px(max(e1a,e1b)) - int(L))*scale, (py(max(e1a,e1b))-int(T))*scale), ((px(min(e1a,e1b))-int(L))*scale, (py(min(e1a,e1b))-int(T))*scale)], fill=(200,0,200), width=1)
crop.save(out)
print(crop.size)

import sys
from PIL import Image
import numpy as np
# usage: axes.py image x0 y0 x1 y1  -> report long dark horizontal/vertical lines in crop (frame candidates)
im=np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int)
x0,y0,x1,y1=map(int,sys.argv[2:6])
c=im[y0:y1,x0:x1]
dark=(c.sum(axis=2)<200)
rows=dark.sum(axis=1); cols=dark.sum(axis=0)
W=x1-x0; H=y1-y0
print("rows with >60% dark:",[ (i+y0,int(r)) for i,r in enumerate(rows) if r>0.6*W])
print("cols with >60% dark:",[ (j+x0,int(r)) for j,r in enumerate(cols) if r>0.6*H])

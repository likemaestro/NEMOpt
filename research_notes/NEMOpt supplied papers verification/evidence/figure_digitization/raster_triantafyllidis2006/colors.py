import sys
from PIL import Image
import numpy as np
from collections import Counter
im=np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int)
x0,y0,x1,y1=map(int,sys.argv[2:6])
c=im[y0:y1,x0:x1].reshape(-1,3)
# quantize and count non-gray colors
cnt=Counter()
for r,g,b in c:
    if max(r,g,b)-min(r,g,b)>60:
        cnt[(r//32*32,g//32*32,b//32*32)]+=1
for k,v in cnt.most_common(15): print(k,v)

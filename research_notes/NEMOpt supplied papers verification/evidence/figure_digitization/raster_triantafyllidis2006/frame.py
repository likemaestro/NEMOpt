import sys
from PIL import Image
import numpy as np
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
H,W,_ = im.shape
dark = (im.sum(axis=2) < 250)
# long horizontal lines
rows = dark.sum(axis=1)
cols = dark.sum(axis=0)
print('size',W,H)
print('rows with >400 dark px:', [ (i,int(rows[i])) for i in range(H) if rows[i]>400])
print('cols with >400 dark px:', [ (i,int(cols[i])) for i in range(W) if cols[i]>400])

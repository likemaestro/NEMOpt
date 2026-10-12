import sys, json
from PIL import Image
import numpy as np
# args: image xL xR yT yB e1L e1R e2T e2B excl_json colors_json angles_deg_json
img=sys.argv[1]
xL,xR,yT,yB=map(float,sys.argv[2:6])
e1L,e1R,e2T,e2B=map(float,sys.argv[6:10])
excl=json.loads(sys.argv[10]); want=json.loads(sys.argv[11]); angs=json.loads(sys.argv[12])
im=np.asarray(Image.open(img).convert("RGB")).astype(int)
Y0,Y1,X0,X1=int(yT)+5,int(yB)-4,int(xL)+5,int(xR)-4
ys,xs=np.mgrid[Y0:Y1,X0:X1]; sub=im[Y0:Y1,X0:X1]
r,g,b=sub[...,0],sub[...,1],sub[...,2]
cls={
 "red":  (r>170)&(g<90)&(b<70)&(r-g>110),
 "green":(g>90)&(r<110)&(g-r>50)&(g-b>30),
 "blue": (b>100)&(r<120)&(g<120)&(b-r>50)&(b-g>50),
 "magenta":(r>170)&(b>90)&(g<100),
 "black":(r+g+b<200),
}
e1=e1L+(xs-xL)/(xR-xL)*(e1R-e1L); e2=e2T+(ys-yT)/(yB-yT)*(e2B-e2T)
for bx in excl:
    m=(e1>=bx[0])&(e1<=bx[1])&(e2>=bx[2])&(e2<=bx[3])
    for k in cls: cls[k]=cls[k]&~m
phi=np.degrees(np.arctan2(e2,e1))%360; rad=np.hypot(e1,e2)
for name in want:
    m=cls[name]; P=phi[m]; R=rad[m]
    row=[]
    for a in angs:
        sel=np.abs(P-a)<0.5
        if sel.sum()==0: row.append((a,None)); continue
        rs=np.sort(R[sel])
        # cluster radii
        cl=[[rs[0]]]
        for v in rs[1:]:
            if v-cl[-1][-1]<0.004: cl[-1].append(v)
            else: cl.append([v])
        row.append((a,[round(float(np.mean(c)),4) for c in cl if len(c)>=2]))
    print(name)
    for a,v in row: print("  phi=%5.1f deg  radius(s)=%s"%(a,v))

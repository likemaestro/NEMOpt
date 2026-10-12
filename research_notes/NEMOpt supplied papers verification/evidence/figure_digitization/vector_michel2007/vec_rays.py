import sys, json, math, pdfplumber
# args: pdf page frame_json(x0,x1,ytop,ybot) limits_json(e1L,e1R,e2T,e2B) angles_json
pdf=pdfplumber.open(sys.argv[1]); pg=pdf.pages[int(sys.argv[2])-1]
fx0,fx1,fyt,fyb=json.loads(sys.argv[3]); e1L,e1R,e2T,e2B=json.loads(sys.argv[4]); angs=json.loads(sys.argv[5])
def to_data(x,y):
    return (e1L+(x-fx0)/(fx1-fx0)*(e1R-e1L), e2T+(y-fyt)/(fyb-fyt)*(e2B-e2T))
def dclass(d):
    if d is None: return "solid"
    arr=d[0] if isinstance(d,(tuple,list)) else []
    if len(arr)==0: return "solid"
    if len(arr)==2: return "dot%.1f/%.1f"%tuple(arr) 
    return "pat"+"/".join("%.1f"%a for a in arr)
groups={}
for o in pg.curves+pg.lines:
    if o["x1"]<fx0-1 or o["x0"]>fx1+1 or o["bottom"]<fyt-1 or o["top"]>fyb+1: continue
    col=o.get("stroking_color"); 
    if col is None: continue
    col=tuple(round(c,2) for c in col) if isinstance(col,(tuple,list)) else col
    if col==(0.0,0.0,0.0) and max(o["x1"]-o["x0"],o["bottom"]-o["top"])>100: continue
    key=(str(col),dclass(o.get("dash")))
    raw=[]
    path=o.get("path") or [("m",q) for q in o["pts"]]
    cur=None
    for cmd in path:
        if cmd[0] in ("m","l"):
            cur=cmd[1]; raw.append(cur)
        elif cmd[0]=="c" and cur is not None:
            p0=cur; p1,p2,p3=cmd[1],cmd[2],cmd[3]
            for k in range(1,21):
                t=k/20; mt=1-t
                raw.append((mt**3*p0[0]+3*mt*mt*t*p1[0]+3*mt*t*t*p2[0]+t**3*p3[0], mt**3*p0[1]+3*mt*mt*t*p1[1]+3*mt*t*t*p2[1]+t**3*p3[1]))
            cur=p3
        elif cmd[0]=="h" and raw:
            raw.append(raw[0])
    pts=[to_data(x,y) for (x,y) in raw]
    groups.setdefault(key,[]).append(pts)
def ray_hits(polys,phi):
    c,s=math.cos(phi),math.sin(phi); hits=[]
    for P in polys:
        for (a1,a2),(b1,b2) in zip(P[:-1],P[1:]):
            # solve t*(c,s) = A + u*(B-A)
            d1,d2=b1-a1,b2-a2
            den=c*d2-s*d1
            if abs(den)<1e-14: continue
            u=(a1*s-a2*c)/den  # from cross products
            t=(a1*d2-a2*d1)/den
            if -1e-9<=u<=1+1e-9 and t>0: hits.append(round(t,4))
    return sorted(set(hits))
for key,polys in sorted(groups.items()):
    npts=sum(len(p) for p in polys)
    if npts<6: continue
    print(key,"paths=%d pts=%d"%(len(polys),npts))
    for a in angs:
        h=ray_hits(polys,math.radians(a))
        if h: print("   phi=%6.1f  lambda=%s   (e1,e2)=%s"%(a,h,[(round(t*math.cos(math.radians(a)),4),round(t*math.sin(math.radians(a)),4)) for t in h[:3]]))

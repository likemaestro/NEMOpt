import sys, pdfplumber
from collections import Counter
pdf=pdfplumber.open(sys.argv[1]); pg=pdf.pages[int(sys.argv[2])-1]
print("page size",pg.width,pg.height)
for kind in ("lines","curves","rects"):
    objs=getattr(pg,kind)
    c=Counter()
    for o in objs:
        c[(str(o.get("stroking_color")),str(o.get("dash")),round(o.get("linewidth",0),2))]+=1
    print(kind,len(objs))
    for k,v in c.most_common(30): print("   ",v,k)
print("images",[(round(i["x0"]),round(i["top"]),round(i["x1"]),round(i["bottom"])) for i in pg.images])
o=(pg.curves or pg.lines)[0]; print(sorted(o.keys()))

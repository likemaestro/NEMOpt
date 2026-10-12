import sys, json, pdfplumber
pdf=pdfplumber.open(sys.argv[1]); pg=pdf.pages[int(sys.argv[2])-1]
fx0,fx1,fyt,fyb=json.loads(sys.argv[3]); e1L,e1R,e2T,e2B=json.loads(sys.argv[4])
want_col=[round(v,2) for v in json.loads(sys.argv[5])]; want_dash=sys.argv[6]
def to_data(x,y): return (round(e1L+(x-fx0)/(fx1-fx0)*(e1R-e1L),4), round(e2T+(y-fyt)/(fyb-fyt)*(e2B-e2T),4))
for o in pg.curves+pg.lines:
    if o["x1"]<fx0-1 or o["x0"]>fx1+1 or o["bottom"]<fyt-1 or o["top"]>fyb+1: continue
    col=o.get("stroking_color")
    if not col or [round(c,2) for c in col]!=want_col: continue
    d=o.get("dash"); arr=d[0] if d else []
    ds="solid" if not arr else "/".join("%.1f"%a for a in arr)
    if want_dash!="any" and ds!=want_dash: continue
    print(ds, [to_data(x,y) for x,y in o["pts"]])

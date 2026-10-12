import sys, pdfplumber
pdf=pdfplumber.open(sys.argv[1]); pg=pdf.pages[int(sys.argv[2])-1]
ytop=float(sys.argv[3]); ybot=float(sys.argv[4])
print("-- long black solid lines (frame candidates)")
for o in pg.lines+pg.curves:
    if o["top"]<ytop or o["bottom"]>ybot: continue
    if str(o.get("stroking_color"))!="(0.0, 0.0, 0.0)": continue
    L=max(o["x1"]-o["x0"], o["bottom"]-o["top"])
    if L>80: print("  x0=%.2f x1=%.2f top=%.2f bottom=%.2f dash=%s"%(o["x0"],o["x1"],o["top"],o["bottom"],o.get("dash")))
print("-- words")
ws=pg.extract_words()
for w in ws:
    if ytop<=w["top"]<=ybot: print("  %-28s x0=%.1f x1=%.1f top=%.1f"%(w["text"],w["x0"],w["x1"],w["top"]))
print("-- black dashed legend samples")
for o in pg.lines+pg.curves:
    if o["top"]<ytop or o["bottom"]>ybot: continue
    if str(o.get("stroking_color"))=="(0.0, 0.0, 0.0)" and o.get("dash") not in (None,([],0)) :
        print("  ",o.get("dash"),"x0=%.1f x1=%.1f top=%.1f"%(o["x0"],o["x1"],o["top"]))

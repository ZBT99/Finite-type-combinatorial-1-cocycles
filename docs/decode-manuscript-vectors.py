from pathlib import Path
import pdfplumber,math,json,hashlib,re
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('gdf_directory')
parser.add_argument('--output',default='manuscript-gdf.json')
args=parser.parse_args()
root=Path(args.gdf_directory)
records={}
for p in sorted(root.glob('*.pdf')):
 with pdfplumber.open(p) as doc:
  pg=doc.pages[0];center=(pg.width/2,pg.height/2);arrows=[]
  for c in pg.curves:
   if not c['stroke'] or c['fill'] or len(c['pts'])!=2:continue
   # Matplotlib's open cubic path runs from tail to the tip of its arrowhead.
   color=c['stroking_color'];req=1 if isinstance(color,tuple) and color[0]>.5 else 0
   arrows.append((c['pts'][0],c['pts'][1],req))
  ends=[]
  for label,(tail,head,req) in enumerate(arrows,1):
   for pt,d in [(tail,-1),(head,1)]:
    angle=(math.atan2(-(pt[1]-center[1]),pt[0]-center[0])-math.pi/2)%(2*math.pi)
    ends.append((angle,[label,d,req,0 if req else 1]))
  ordered=[r for a,r in sorted(ends)];mapping={}
  for r in ordered:mapping.setdefault(r[0],len(mapping)+1);r[0]=mapping[r[0]]
  records[p.stem]=dict(rows=ordered,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
Path(args.output).write_text(json.dumps(records,indent=2))
print(len(records),'vector PDFs decoded')

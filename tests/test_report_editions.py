"""Language parity and arrowhead readability regressions."""
from html.parser import HTMLParser
from pathlib import Path
import math,re,unittest
import numpy as np
from cocycle.diagrams import diagram
from cocycle.language import LANG
from cocycle.reports import page
ROOT=Path(__file__).resolve().parents[1]

class ReportParser(HTMLParser):
    def __init__(self):
        super().__init__();self.text=[];self.tags=[];self.geometry=[];self.pre=[];self.in_pre=False;self.buffer=[];self.in_td=False;self.cell=[];self.numeric_cells=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs);self.tags.append(tag)
        if tag=='pre':self.in_pre=True;self.buffer=[]
        if tag=='td':self.in_td=True;self.cell=[]
        if 'data-arrow' in d:self.geometry.append(('arrow',d['data-arrow'],d['d'],d['stroke'],d['stroke-width']))
        if 'data-endpoint' in d:self.geometry.append(('dot',d['data-endpoint'],d['cx'],d['cy'],d['r']))
    def handle_endtag(self,tag):
        if tag=='pre':self.pre.append(''.join(self.buffer));self.in_pre=False
        if tag=='td':
            value=''.join(self.cell).strip()
            if re.fullmatch(r'[0-9\s+\-−×()%.,\[\]]+',value):self.numeric_cells.append(value)
            self.in_td=False
    def handle_data(self,data):
        self.text.append(data)
        if self.in_pre:self.buffer.append(data)
        if self.in_td:self.cell.append(data)

class ReportEditionChecks(unittest.TestCase):
    def test_all_editions_have_same_structure_geometry_and_raw_data(self):
        english=sorted((ROOT/'reports/en').glob('*.html'))
        self.assertEqual(len(english),11)
        self.assertEqual({p.name for p in english},{p.name for p in (ROOT/'reports/zh').glob('*.html')})
        for p in english:
            en,zh=ReportParser(),ReportParser();en.feed(p.read_text());zh.feed((ROOT/'reports/zh'/p.name).read_text())
            self.assertIsNone(re.search(r'[\u3400-\u9fff]',''.join(en.text)),p.name)
            self.assertIsNotNone(re.search(r'[\u3400-\u9fff]',''.join(zh.text)),p.name)
            # Translated grammar can reorder inline emphasis, but not report sections or tables.
            structure={'h1','h2','h3','p','details','summary','table','tr','th','td','li','pre','svg'}
            self.assertEqual([t for t in en.tags if t in structure],[t for t in zh.tags if t in structure],p.name)
            self.assertEqual(en.numeric_cells,zh.numeric_cells,p.name)
            self.assertEqual(en.geometry,zh.geometry,p.name)
            self.assertEqual(en.pre,zh.pre,p.name)

    def test_gray_arrowhead_is_fixed_size_and_clear_of_endpoint(self):
        g=np.array([[1,-1,1],[2,-1,-1],[1,1,1],[2,1,-1]])
        markup=diagram(g)
        self.assertIn('markerUnits="userSpaceOnUse"',markup)
        self.assertIn('markerWidth="11" markerHeight="9"',markup)
        parser=ReportParser();parser.feed(markup)
        dots={int(v[1]):(float(v[2]),float(v[3])) for v in parser.geometry if v[0]=='dot'}
        for kind,label,path,stroke,width in parser.geometry:
            if kind!='arrow':continue
            tip=re.search(r'L([\d.-]+),([\d.-]+)',path)
            end=int(np.flatnonzero((g[:,0]==int(label))&(g[:,1]==1))[0]);x,y=dots[end]
            distance=math.hypot(x-float(tip[1]),y-float(tip[2]))
            self.assertAlmostEqual(distance,7,delta=.02)
            self.assertEqual(width,'1.5')
        self.assertIn('stroke="#53606e"',markup)

    def test_language_context_restored(self):
        self.assertEqual(LANG.get(),'en')
        self.assertIn('lang="zh-CN"',page('测试','',lang='zh'))
        self.assertEqual(LANG.get(),'en')
        with self.assertRaises(ValueError):page('','',lang='xx')

if __name__=='__main__':unittest.main()

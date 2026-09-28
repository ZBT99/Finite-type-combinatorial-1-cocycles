"""Read-only localhost explorer with complete English and Chinese reports."""
from http.server import HTTPServer,BaseHTTPRequestHandler
from urllib.parse import urlparse,parse_qs,urlencode
from pathlib import Path
import json
from .reports import page,esc,local_case,commutation_case,loop_case,para,heading
from .language import tr,localized,LANG
from .loops import knots,calculate

@localized
def home(bundle,loopresults):
    options=''.join('<option>'+esc(n)+'</option>' for n in knots())
    beta='<label>Cocycle <select name="beta"><option>beta1</option><option>beta2</option><option>beta3</option></select></label>'
    hidden=f'<input type="hidden" name="lang" value="{LANG.get()}">'
    button='<p><button>'+tr('Show all contributions','查看全部贡献')+'</button></p>'
    body=para('Select a case to inspect all occurrences, signs, coefficients and the exact identity certificate.','选择 case，查看全部匹配、符号、系数及精确恒等式证书。')
    body+=heading('Local equations','局部方程')+para('Tetrahedron and cube equations vanish on classical diagrams by the displayed linking identities. Raw virtual residuals remain visible.','Tetrahedron 与 cube 方程通过展示的 linking identities 在经典图上为零。虚拟图的原始残差仍会显示。')+'<div class="grid">'
    body+='<form class="card" action="/case">'+hidden+'<h3>Tetrahedron</h3><input hidden name="family" value="tetra"><p>s (1–24) <input name="case" type="number" min="1" max="24" value="1"></p>'+beta+button+'</form>'
    body+='<form class="card" action="/case">'+hidden+'<h3>Cube</h3><input hidden name="family" value="cube"><p>'+tr('Template','模板')+' (0–47) <input name="template" type="number" min="0" max="47" value="1"></p><p>'+tr('Basepoint rotation','基点旋转')+' (0–2) <input name="rotation" type="number" min="0" max="2" value="1"></p>'+beta+button+'</form>'
    body+='<form class="card" action="/case">'+hidden+'<h3>R3–R3 commute</h3><input hidden name="family" value="commute"><p>'+tr('First type','第一类型')+' (0–47) <input name="u" type="number" min="0" max="47" value="0"></p><p>'+tr('Second type','第二类型')+' (0–47) <input name="v" type="number" min="0" max="47" value="0"></p><p>'+tr('Interleaving','交错方式')+' (0–9) <input name="a" type="number" min="0" max="9" value="0"></p>'+beta+button+'</form></div>'
    body+=heading('Loops and push arcs','Loops 与 push arcs')+'<form action="/loop">'+hidden+'<p><select name="kind">'+''.join('<option>'+k+'</option>' for k in ['rotation','rolling','half-rolling','bracket','half-bracket','push'])+'</select> <select name="knot">'+options+'</select></p><p>'+tr('Second knot (push/bracket only)','第二个结（仅 push/bracket 使用）')+': <select name="other"><option value=""></option><option>curl+</option><option>2curls+</option>'+options+'</select></p>'
    body+=para('Use 4_1_periodic or 6_1_periodic for half rolling. Rolling follows the paper orientation; the raw notebook direction is shown separately.','Half rolling 使用 4_1_periodic 或 6_1_periodic。Rolling 遵循论文方向，另行显示原 notebook 方向。')+'<button>'+tr('Compute and display','计算并显示')+'</button></form>'
    body+=heading('Verified results','验证结果')+para('Local cores: 24 tetrahedron and 144 cube. R3–R3 cases: {n}. Loop comparisons: {m}.','局部核心：24 个 tetrahedron、144 个 cube。R3–R3 情况：{n}。Loop 比较：{m}。',n=bundle['dr3']['cases'],m=loopresults['cases'])
    body+=para('β₃ minimal combinatorial order is outside this package. Loop examples supply computational input; the general theorems also rely on the paper’s mathematical arguments.',
               '本包不证明 β₃ 的最小组合阶数。Loop 例子提供计算输入；一般定理还依赖论文中的数学论证。')
    return page(tr('Order-four cocycles: case explorer','四阶 cocycle：情况浏览器'),body)

def serve(port,results,lang='en'):
    results=Path(results);bundle=json.loads((results/'local-certificates.json').read_text());loopresults=json.loads((results/'loop-checks.json').read_text())
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            url=urlparse(self.path);q={k:v[0] for k,v in parse_qs(url.query).items()};selected=q.get('lang',lang)
            if selected not in ('en','zh'):self.send_error(400);return
            token=LANG.set(selected)
            try:
                if url.path=='/':content=home(bundle,loopresults)
                elif url.path=='/guide':
                    from .examples import reading_guide
                    content=reading_guide(bundle).replace('cube-1-1-beta1.html#case-16-18',f'/case?family=cube&amp;template=1&amp;rotation=1&amp;beta=beta1&amp;lang={selected}#case-16-18')
                elif url.path=='/case':
                    family=q['family'];name=q.get('beta','beta1')
                    if name not in ['beta1','beta2','beta3']:raise ValueError('Unknown cocycle')
                    if family=='commute':
                        u,v,a=(int(q[k]) for k in ('u','v','a'))
                        if not (0<=u<48 and 0<=v<48 and 0<=a<10):raise ValueError('Case index outside range')
                        content=commutation_case(u,v,a,name=name)
                    else:content=local_case(bundle,family,int(q['case']) if family=='tetra' else (int(q['template']),int(q['rotation'])),name)
                elif url.path=='/loop':content=loop_case(calculate(q['kind'],q['knot'],q.get('other') or None))
                else:self.send_error(404);return
            except (ValueError,KeyError,IndexError,StopIteration):
                content=page(tr('Invalid case','无效的 case'),para('Check the case indices and required knot inputs, then return to the selector.','请检查 case 索引和必需的结输入，再返回选择页面。'))
            finally:LANG.reset(token)
            content=content.replace('href="reading-guide.html"',f'href="/guide?lang={selected}"').replace('href="index.html"',f'href="/?lang={selected}"')
            other='zh' if selected=='en' else 'en';alternate=dict(q,lang=other)
            link=url.path+'?'+urlencode(alternate);label='Chinese edition' if selected=='en' else '英文版'
            content=content.replace('</nav>',f'<a href="{esc(link)}">{label}</a></nav>',1)
            self.send_response(200);self.send_header('Content-Type','text/html; charset=utf-8');self.end_headers();self.wfile.write(content.encode())
    print(f'Open http://127.0.0.1:{port}/?lang={lang} (Ctrl-C to stop)',flush=True)
    HTTPServer(('127.0.0.1',port),Handler).serve_forever()

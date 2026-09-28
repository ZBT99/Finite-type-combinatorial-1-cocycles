"""Build complete, parallel English and Chinese editions of every offline report."""
import base64,json
from pathlib import Path
import numpy as np
from .reports import page,local_case,commutation_case,loop_case,para,heading,esc,beta_label
from .diagrams import diagram,arc_spec
from .language import tr,localized,LANG
from .verification import dr3_case
from . import local as b,invariants

@localized
def reading_guide(bundle):
    cards=[]
    samples=[('Cube',bundle['cube'][0]['core'],bundle['cube'][0]['indices'][0],arc_spec('cube',0)),
             ('Tetrahedron',bundle['tetra'][0]['core'],bundle['tetra'][0]['indices'][0],arc_spec('tetra'))]
    gs,ix=dr3_case(0,2,3);samples.append(('R3–R3 commute',gs[0],ix[0],arc_spec('commute')))
    for title,g,index,spec in samples:
        m=len(spec['sizes'])
        cards.append('<section><h3>'+title+tr(' · {m} solid blocks',' · {m} 段实线块',m=m)+'</h3>'+diagram(g,index,spec=spec)+para('Endpoints per core block: {s}. Dashed gaps: B<sub>0</sub>–B<sub>{m}</sub>.','每个核心块的端点数：{s}。虚线补区间：B<sub>0</sub>–B<sub>{m}</sub>。',s=list(spec['sizes']),m=m)+'</section>')
    body=para('Read the arc labels first, then arrow directions and crossing signs, and finally the residual and its certificate.','先看圆弧编号，再看箭头方向和交叉符号，最后看残差及其证书。')
    body+=heading('1. The three core partitions','1. 三种核心分块')
    body+=para('Solid circle arcs connect endpoints in a fixed core block. Dashed arcs are gaps where outside endpoints may be inserted. Filled dots are core endpoints; hollow dots are outside. Gaps stay dashed after insertion. The partition belongs to the whole core, not just the current R3 triangle.',
               '实线圆弧连接固定核心块中的端点；虚线圆弧是可插入外部端点的补区间。实心点是核心端点，空心点是外部端点；插入后补弧仍为虚线。这套分块属于整个核心，不仅属于当前 R3 三角形。')
    body+=para('Read counterclockwise from ∞ at the top: B<sub>0</sub> lies to its left, and the last B lies to its right. The basepoint splits one circular gap into two based intervals, so m core blocks give m+1 gap labels.',
               '从顶部 ∞ 逆时针读取：左侧是 B<sub>0</sub>，右侧是最后一个 B。基点将一个圆周间隙分成两个 based intervals，因此 m 个核心块对应 m+1 个补区间编号。')
    body+='<div class="grid">'+''.join(cards)+'</div>'
    body+=para('Commute labels B<sub>0</sub>–B<sub>6</sub> are only a reading aid; no additional commute outside-arrow enumeration is claimed. The four intervals in an individual R3 formula array are a different convention from these whole-core B labels.',
               'Commute 的 B<sub>0</sub>–B<sub>6</sub> 只辅助读图，不表示程序额外枚举了外部箭头。单个 R3 公式数组里的四个区间，与这里整个核心的 B 编号是不同约定。')
    body+=heading('2. Arrowheads and endpoints','2. 箭头头部与端点')
    body+=para('Arrows point from underpass to overpass. Their arrowheads have a fixed display size, independent of shaft thickness, and stop before the endpoint dot. Red marks the current R3; blue marks selected contributing arrows; gray marks the remaining arrows. A label such as 5:− gives crossing number and sign. Intersecting chords are not new knot crossings.',
               '箭头从下穿端指向上穿端。箭头头部大小不依赖线条粗细，并停在端点圆点之前。红色标当前 R3，蓝色标选中贡献箭头，灰色标其余箭头。标签 5:− 表示交叉编号和符号。弦相交不是新的结交叉。')
    body+=heading('3. What X[a,b,d] counts','3. X[a,b,d] 统计什么')
    body+=para('a,b are gap numbers with a ≤ b, not crossing labels. At the endpoint encountered first from ∞, d=+1 means head/overpass and d=−1 means tail/underpass. When a=b, use the order within that gap. <b>d is not the crossing sign.</b> Each matching positive crossing adds +1 to X; each negative crossing adds −1.',
               'a、b 是补区间编号，a ≤ b，不是交叉编号。从 ∞ 出发先遇到的端点若为箭头端／上穿端则 d=+1，若为箭尾端／下穿端则 d=−1。a=b 时按同一补区间内的顺序判断。<b>d 不是交叉符号。</b>满足条件的正交叉向 X 加 +1，负交叉加 −1。')
    body+=para('Thus X[1,2,−1] sums signs of arrows from B<sub>1</sub> to B<sub>2</sub>. In X[1,2,−1]=−1, the −1 inside brackets specifies direction; the −1 after the equals sign is the signed count. The subscript b in B<sub>b</sub> is a Latin letter. The constant term is a signed core contribution, not a count of core arrows.',
               '因此 X[1,2,−1] 是从 B<sub>1</sub> 指向 B<sub>2</sub> 的箭头的交叉符号之和。在 X[1,2,−1]=−1 中，括号里的 −1 表方向，等号后的 −1 是带符号计数。B<sub>b</sub> 的下标 b 是拉丁字母。常数项是核心交叉的带符号贡献，不是核心箭头根数。')
    body+=heading('4. Residual increment','4. 残差增量')
    body+='<div class="equation">'+tr('ΔR = R(core + one outside arrow) − R(core)','ΔR = R(核心＋一根外部箭头) − R(核心)')+'</div>'
    body+=para('The heading shows only the selected cocycle. Core residual plus increment equals the full one-arrow residual. β₃ is reduced modulo two, with its integer lift retained. Event values already include coorientation and path sign; do not apply these factors again.',
               '标题只显示所选 cocycle。核心残差加增量，等于完整单箭头图的残差。β₃ 按模二计算，并保留整数提升。事件值已乘 coorientation 和路径方向，不要重复乘这些因子。')
    body+=heading('5. Cube cases 16 and 18','5. Cube 的 case 16 与 18')
    body+=para('For cube (1,1), case 16 has a negative arrow B<sub>1</sub> → B<sub>2</sub>: term 6 at p2 gives −1. Case 18 has a negative arrow in the opposite direction: term 10 at p1 gives +1. Here the core residual is zero.',
               '对 cube (1,1)，case 16 有一根负交叉箭头 B<sub>1</sub> → B<sub>2</sub>：p2 的第 6 项给出 −1。Case 18 的负交叉箭头方向相反：p1 的第 10 项给出 +1。这里核心残差为零。')
    body+='<div class="equation">R = X[1,2,−1] − X[1,2,+1]</div>'
    body+=para('One arrow of each type in a single complete diagram would cancel these increments. This is not a proof by summing independent test cases, and it does not construct a classical completion. Smoothing crossings 2,3 makes crossing 5 the only crossing between components 1,2 in either isolated one-arrow case; this cannot be a classical two-component diagram.',
               '若同一张完整图中两种箭头各有一根，它们的增量相消。这不是把独立测试 case 相加作为证明，也没有构造出一张经典补全图。平滑交叉 2、3 后，这两个孤立的单箭头 case 中，交叉 5 都是分量 1、2 之间唯一的交叉；这不可能是经典两分量图。')
    body+=heading('6. How the certificate proves classical cancellation','6. 证书如何证明经典情形的整体消去')
    steps=[('Compute R from actual R3 sums, including the constant and every outside-arrow coefficient.','实际计算 R3 求和，得到 R 的常数项与每个外部箭头系数。'),
           ('Independently smooth core crossings, trace components, and construct linking identities H from signed over/under counts.','独立平滑核心交叉、追踪分量，根据上下穿的带符号计数构造 linking identities H。'),
           ('Check every coefficient exactly to verify R = Σ cᵢHᵢ.','精确逐系数核对 R = Σ cᵢHᵢ。')]
    body+='<ol>'+''.join('<li>'+tr(*s)+'</li>' for s in steps)+'</ol>'
    body+=para('For classical links, Hᵢ=0 because the over- and under-crossing sums calculate the same linking number. Hence R=0 for classical completions. Formal virtual cases can have nonzero raw residuals. No unique one-to-one arrow matching is asserted. Arithmetic is exact over Q for β₁,β₂ and over F₂ for β₃.',
               '对经典链环，上下穿的两个带符号交叉和计算同一个 linking number，因此 Hᵢ=0，从而经典补全图满足 R=0。形式虚拟图可以有非零原始残差。这里不声称存在唯一的逐根箭头配对。β₁、β₂ 在 Q 上精确运算，β₃ 在 F₂ 上运算。')
    body+='<p><a href="cube-1-1-beta1.html#case-16-18">'+tr('Open the worked cube example','打开 cube 具体图例')+'</a></p>'
    return page(tr('Reading diagrams and certificates','读图与证书指南'),body)

@localized
def formula_catalogue(root):
    body=para('The β figures reproduce the supplied manuscript images. Term numbering follows the arrays. A configuration row is (arrow label, endpoint direction, required crossing sign, R3 interval). A zero required sign is unrestricted. Unsigned terms multiply crossing signs; signed terms require their displayed signs and have weight one.',
              'β 公式图使用提供的论文原图。项编号遵循数组顺序。配置每行为（箭头编号、端点方向、要求的交叉符号、R3 区间）。符号要求为 0 表示不限制。未指定符号的项乘交叉符号；指定符号的项要求相应符号，权重为 1。')
    for i,name in enumerate(b.bv_names,1):
        f=b.bv_formulas[name];image=base64.b64encode((root/f'docs/formulas/beta_{i}.png').read_bytes()).decode()
        body+=heading(beta_label(name),beta_label(name))+f'<img style="max-width:100%" alt="{beta_label(name)}" src="data:image/png;base64,{image}">'
        body+='<div class="table-wrap"><table><tr>'+''.join('<th>'+tr(*h)+'</th>' for h in [('Term','项'),('Coefficient','系数'),('R3 type','R3 类型'),('Configuration rows','配置各行')])+'</tr>'
        for j,(d,c,t) in enumerate(zip(f['cfg'],f['coeff'],f['types']),1):
            rows=d[d[:,0]!=0].tolist();body+=f'<tr><td>{j}</td><td>{int(c)}</td><td>{int(t)}</td><td><code>{esc(rows)}</code></td></tr>'
        body+='</table></div>'
    body+=heading('Invariant formulas','不变量公式')
    body+=para('For these invariants, the fourth column is a sign exponent, not an R3 interval. An endpoint row is (label, direction, required sign, exponent). A chord marked ± has unrestricted crossing sign; consult the exponent in the rows for its weight. v₄,₁=J15 and v₄,₂=E34.',
               '这些不变量公式的第四列是交叉符号的指数，不是 R3 区间。每行为（编号、方向、要求的符号、指数）。弦上 ± 表示交叉符号不受限制，具体权重由各行的指数决定。v₄,₁=J15，v₄,₂=E34。')
    for name,label in [('v2','v₂'),('v3','v₃'),('v41','v₄,₁ = J15'),('v42','v₄,₂ = E34')]:
        ds,cs=invariants.formula(name);body+=heading(label,label)+'<div class="grid">'
        for j,(d,c) in enumerate(zip(ds,cs),1):
            d=d[d[:,0]!=0]
            body+='<div class="card"><h3>'+tr('Term {j}; coefficient {c}','第 {j} 项；系数 {c}',j=j,c=int(c))+'</h3>'+diagram(d[:,:3])+ '<pre>'+esc(json.dumps(d.tolist()))+'</pre></div>'
        body+='</div>'
    return page(tr('Formula catalogue','公式目录'),body)

@localized
def example_index():
    body=para('All reports work offline without installation or remote scripts. Start with the reading guide, then inspect the worked cube example. The two language editions contain the same calculations and corresponding explanatory paragraphs.',
              '所有报告均可离线打开，无需安装或远程脚本。建议先读指南，再看 cube 图例。两种语言版本使用相同计算数据，解释段落逐一对应。')
    links=[('reading-guide.html','Reading guide','读图指南'),('cube-1-1-beta1.html#case-16-18','Cube cases 16 and 18 explained','Cube 的 case 16 与 18 详解'),('formulas.html','Formula catalogue','公式目录'),
           ('tetra-1-beta1.html','Tetrahedron 1 · β₁','Tetrahedron 1 · β₁'),('tetra-2-beta2.html','Tetrahedron 2 · β₂','Tetrahedron 2 · β₂'),('tetra-17-beta3.html','Tetrahedron 17 · β₃','Tetrahedron 17 · β₃'),('cube-1-1-beta1.html','Cube (1,1) · β₁ · all outside cases','Cube (1,1) · β₁ · 全部外部箭头 case'),('commute-beta2.html','R3–R3 commute (0,2,3) · β₂','R3–R3 commute (0,2,3) · β₂'),('rolling-4_1.html','Rolling · 4₁','Rolling · 4₁'),('half-rolling-6_1.html','Half rolling · periodic 6₁','Half rolling · 周期 6₁'),('bracket-4_1-3_1.html','Bracket · 4₁, 3₁⁺','Bracket · 4₁、3₁⁺')]
    body+='<ul>'+''.join(f'<li><a href="{href}">'+tr(en,zh)+'</a></li>' for href,en,zh in links)+'</ul>'
    body+=heading('Scope','验证范围')
    body+=para('Exact classical certificates cover 24 positive tetrahedron cores and 144 cube cores. All 23,040 signed R3–R3 meridian sums vanish; 136 loop comparisons agree. Raw virtual residuals are retained. β₃ minimal-order proof is excluded; general loop formulas also use the theoretical arguments in the paper.',
               '精确经典证书覆盖 24 个正 tetrahedron 核心和 144 个 cube 核心。23,040 个带符号 R3–R3 meridian 总和均为零；136 次 loop 比较一致。虚拟图的原始非零残差被保留。本包不证明 β₃ 的最小阶数；一般 loop 公式还依赖论文中的理论论证。')
    body+=para('To select other cases, run <code>python -m cocycle serve --lang en</code> and open <code>http://127.0.0.1:8765</code>. Rebuild both editions with <code>python -m cocycle.examples</code>.',
               '选择其他情况时，运行 <code>python -m cocycle serve --lang zh</code>，打开 <code>http://127.0.0.1:8765</code>。用 <code>python -m cocycle.examples</code> 重新生成两种语言版本。')
    return page(tr('Order-four cocycles: example reports','四阶 cocycle：示例报告'),body)


def edition_link(document,filename,lang):
    other='zh' if lang=='en' else 'en'
    label='Chinese edition' if lang=='en' else '英文版'
    return document.replace('</nav>',f'<a href="../{other}/{filename}">{label}</a></nav>',1)


def main():
    root=Path.cwd();out=root/'reports';bundle=json.loads((root/'results/local-certificates.json').read_text())
    entries=[('tetra-1-beta1.html','tetra',1,'beta1'),('tetra-2-beta2.html','tetra',2,'beta2'),('tetra-17-beta3.html','tetra',17,'beta3'),('cube-1-1-beta1.html','cube',(1,1),'beta1')]
    for lang in ('en','zh'):
        target=out/lang;target.mkdir(exist_ok=True)
        def save(filename,content):
            (target/filename).write_text(edition_link(content,filename,lang),encoding='utf-8')
            print(lang,filename,flush=True)
        for filename,family,case,name in entries:save(filename,local_case(bundle,family,case,name,lang=lang))
        save('commute-beta2.html',commutation_case(0,2,3,name='beta2',lang=lang))
        for stem in ('rolling-4_1','half-rolling-6_1','bracket-4_1-3_1'):
            save(stem+'.html',loop_case(json.loads((out/(stem+'.json')).read_text()),lang=lang))
        save('reading-guide.html',reading_guide(bundle,lang=lang))
        save('formulas.html',formula_catalogue(root,lang=lang))
        save('index.html',example_index(lang=lang))
    # Keep the root as a language chooser; actual reports have one language each.
    for old in out.glob('*.html'):
        if old.name!='index.html':old.unlink()
    chooser='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Report language / 报告语言</title><style>body{font:20px/1.7 system-ui;max-width:680px;margin:60px auto;padding:24px}a{display:block;margin:24px 0;color:#176c9c}</style><body><h1>Report language / 报告语言</h1><a href="en/index.html">English reports</a><a href="zh/index.html">中文版报告</a><p>Version 1.2.0 · Same computations in both editions.<br>两种语言版本使用相同的计算数据。</p></body></html>'''
    (out/'index.html').write_text(chooser,encoding='utf-8')

if __name__=='__main__':main()

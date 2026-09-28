"""Single-language reports with parallel English and Chinese content."""
import html,json
import numpy as np
from . import local as b,matching
from .verification import dr3_case
from .diagrams import diagram,arc_spec,outside_descriptor,COLORS
from .language import tr,localized,LANG
from .report_style import STYLE

def esc(x):return html.escape(str(x))
def beta_label(name):return {'beta1':'β₁','beta2':'β₂','beta3':'β₃','alpha3':'α₃¹'}.get(name,name)
def para(en,zh,**values):return '<p>'+tr(en,zh,**values)+'</p>'
def heading(en,zh,level=2,id=None):return f'<h{level}'+(f' id="{id}"' if id else '')+'>'+tr(en,zh)+f'</h{level}>'
def value_text(value,name):
    v=int(value)
    return tr('{r} (mod 2; integer lift {v:+d})','{r}（模 2；整数提升 {v:+d}）',r=v%2,v=v) if name=='beta3' else str(v)

def arc_guide(spec,commute=False):
    m=len(spec['sizes'])
    out='<div class="reading-guide">'+heading('How to read the circle','圆弧怎么读')
    out+=para('<b>{m} solid core blocks</b>; {n} dashed based gaps B<sub>0</sub>–B<sub>{m}</sub>. Core endpoints per block: {sizes}.',
              '<b>{m} 段实线核心弧</b>；{n} 段虚线补区间 B<sub>0</sub>–B<sub>{m}</sub>。每个核心块的端点数：{sizes}。',m=m,n=m+1,sizes=esc(list(spec['sizes'])))
    out+=para('Solid arcs connect endpoints within the same fixed core block. Dashed arcs are gaps where outside endpoints may be inserted. Core endpoints are filled dots; outside endpoints are hollow. A gap remains dashed after insertion.',
              '实线弧连接同一个固定核心端点块；虚线弧是可插入外部端点的补区间。核心端点为实心点，外部端点为空心点；插入后补弧仍画虚线。')
    out+=para('Read <b>counterclockwise</b> from ∞ at the top. B<sub>0</sub> is immediately to its left; B<sub>{m}</sub> is immediately to its right before returning. The basepoint splits one circular gap into two based intervals. B subscripts are gap numbers, not crossing labels or the Greek letter β.',
              '从顶部 ∞ <b>逆时针</b>读取：左侧首先经过 B<sub>0</sub>，右侧返回前经过 B<sub>{m}</sub>。基点将一个圆周间隙分成两个 based intervals。B 的下标是补区间编号，不是交叉编号，也不是希腊字母 β。',m=m)
    out+=para('Red chords belong to the current R3; blue chords are selected contributing arrows. Arrows point from underpass to overpass. Arrowheads stop just before endpoint dots so that direction stays visible. A chord label such as 5:− means crossing 5 has negative sign; small outer numbers are zero-based endpoint row indices. Intersecting chords do not represent additional crossings.',
              '红色弦属于当前 R3；蓝色弦是选中的贡献箭头。箭头从下穿端指向上穿端；箭头头部停在端点圆点之前，避免被遮挡。弦标签 5:− 表示交叉 5 的符号为负；圆外小数字是从 0 开始的端点行号。弦相交不代表额外的结交叉。')
    if commute:out+=para('The commute B labels only aid reading; no extra-arrow enumeration is asserted here.','Commute 的 B 编号仅用于读图；这里没有额外枚举外部箭头。')
    return out+'</div>'

@localized
def page(title,body):
    nav='<nav class="jump-links"><a href="reading-guide.html">'+tr('Reading guide','读图指南')+'</a><a href="index.html">'+tr('Examples','示例目录')+'</a></nav>'
    controls='<button onclick="document.querySelectorAll(\'details\').forEach(x=>x.open=true)">'+tr('Expand all','展开全部')+'</button> <button onclick="document.querySelectorAll(\'details\').forEach(x=>x.open=false)">'+tr('Collapse all','收起全部')+'</button>'
    script="""<script>function revealHash(){const el=document.getElementById(decodeURIComponent(location.hash.slice(1)));if(!el)return;let p=el;while(p){if(p.tagName==='DETAILS')p.open=true;p=p.parentElement;}el.scrollIntoView();}addEventListener('hashchange',revealHash);addEventListener('load',revealHash);</script>"""
    return '<!doctype html><html lang="'+('zh-CN' if LANG.get()=='zh' else 'en')+'"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>'+esc(title)+'</title><style>'+STYLE+'</style><body><div class="tag">'+tr('ORDER FOUR · EXACT VERIFICATION','四阶 cocycle · 精确验证')+' · 1.2.0</div>'+nav+'<h1>'+esc(title)+'</h1>'+controls+body+script+'</body></html>'

@localized
def event_report(g,index,scale=1,name=None,number=1,cfgs=None,types=None,coeff=None,names=None,spec=None):
    if cfgs is None:cfgs,types,coeff,names=b.bv_cfgs,b.bv_types,b.bv_coeff_matrix,b.bv_names
    terms,occ=matching.event_terms(np.array(g),np.array(index),cfgs,types,True)
    values=terms@coeff*scale
    core=np.array(g)[np.column_stack((index,np.array(index)+1)).ravel()]
    typ=int(b.bv_basic.r_l_R3_type(core));co=int(b.bv_basic.get_R3_sign(core))
    cards=[];rows=[]
    for o in occ:
        glob=o['term']-1;k=int(np.flatnonzero(coeff[glob])[0]);nm=names[k]
        if name and name!=nm:continue
        start=int(np.flatnonzero(coeff[:,k])[0]);term=glob-start+1;c=int(coeff[glob,k]);value=scale*c*o['value']
        label=tr('{nm}, term {term}','{nm}，第 {term} 项',nm=beta_label(nm),term=term)
        rows.append(f'<tr><td>{label}</td><td>{esc(o["arrows"])}</td><td>{esc(o["crossing_signs"])}</td><td>{c} × {co} × {o["matching_weight"]} × {scale}</td><td>{value}</td></tr>')
        note=tr('Required crossing signs; weight 1','要求指定交叉符号；权重为 1') if o['signed'] else tr('Unsigned arrows; multiply crossing signs','未指定符号的箭头；权重为交叉符号之积')
        cards.append('<div class="card">'+diagram(g,index,o['arrows'],spec=spec)+f'<b>{label}: {value:+d}</b>'+para('Selected arrows: {arrows}','选中箭头：{arrows}',arrows=esc(o['arrows']))+'<span class="muted">'+note+'</span></div>')
    val='; '.join(beta_label(n)+' = '+value_text(v,n) for n,v in zip(names,values) if not name or name==n)
    heads=[tr('Formula term','公式项'),tr('Embedding','选中箭头'),tr('Crossing signs','交叉符号'),tr('coefficient × coorientation × weight × path sign','系数 × coorientation × 权重 × 路径方向'),tr('Contribution','贡献')]
    table='<table><tr>'+''.join('<th>'+h+'</th>' for h in heads)+'</tr>'+''.join(rows)+'</table>' if rows else para('No matching occurrence for the selected formula. All contributions are zero.','所选公式没有匹配项；所有贡献均为零。')
    rawterms={n:(terms[coeff[:,k]!=0]*coeff[coeff[:,k]!=0,k]*scale).tolist() for k,n in enumerate(names) if not name or n==name}
    summary=tr('p{number}: {val} · {count} occurrences · R3 type {typ} · coorientation {co:+d} · path sign {scale:+d}',
               'p{number}：{val} · {count} 个匹配 · R3 类型 {typ} · coorientation {co:+d} · 路径方向 {scale:+d}',number=number,val=val,count=len(cards),typ=typ,co=co,scale=scale)
    return '<details><summary>'+summary+'</summary>'+diagram(g,index,spec=spec)+para('R3 endpoint indices (zero based): {indices}. Read counterclockwise from ∞. Coorientation and path sign are separate factors; displayed contributions already include both. Do not reverse p2 a second time.',
        'R3 端点索引（从 0 开始）：{indices}。从 ∞ 逆时针读取。Coorientation 与路径方向是不同因子；显示贡献已包含二者，不要再对 p2 取负。',indices=esc(list(index)))+'<div class="table-wrap">'+table+'</div><div class="grid">'+''.join(cards)+'</div><details><summary>'+tr('Complete Gauss matrix and all term totals, including zeros','完整 Gauss 矩阵与各项总和（包括零项）')+'</summary><pre>'+esc(json.dumps(dict(gauss=np.array(g).tolist(),term_totals=rawterms),indent=2))+'</pre></details></details>'

@localized
def certificate_report(r,name,variables,spec=None):
    k=b.bv_names.index(name);mod=2 if name=='beta3' else None
    poly=lambda row:b.bv_polynomial_text(row,mod,variables)
    target=np.array(r['target'])[:,k];weights=r['certificates'][name]
    out=[heading('Classical identity certificate','经典情形的整体证书',id='certificate')]
    out.append(para('<b>Individual one-arrow cases need not have zero residual. Independent test cases are not added together.</b> This certificate identifies the residual polynomial with a combination of classical linking identities.',
                    '<b>单箭头 case 的残差不一定为零；程序也不会把独立测试 case 相加。</b>证书验证残差多项式等于经典 linking identities 的线性组合。'))
    out.append('<details open><summary>'+tr('What does X[a,b,d] mean?','X[a,b,d] 是什么意思？')+'</summary>')
    out.append(para('X[a,b,d] is the <b>sum of crossing signs</b> of outside arrows with endpoints in B<sub>a</sub>,B<sub>b</sub> (a ≤ b). Direction d is measured at the endpoint encountered first from ∞: +1 is a head/overpass, −1 a tail/underpass. <b>d is not the crossing sign.</b> For a=b, use the order within that gap.',
                    'X[a,b,d] 是端点位于 B<sub>a</sub>、B<sub>b</sub>（a ≤ b）的外部箭头的<b>交叉符号之和</b>。d 表示从 ∞ 出发先遇到的端点方向：+1 是箭头端／上穿端，−1 是箭尾端／下穿端。<b>d 不是交叉符号。</b>若 a=b，则按同一补区间内的端点顺序判断。'))
    out.append(para('A matching positive crossing adds +1 to X; a negative crossing adds −1. For example, X[1,2,−1] counts arrows from B<sub>1</sub> to B<sub>2</sub> with these signs. The constant records the signed contribution of core crossings, not their number.',
                    '满足条件的正交叉向 X 贡献 +1，负交叉贡献 −1。例如 X[1,2,−1] 统计从 B<sub>1</sub> 指向 B<sub>2</sub> 的箭头的带符号总和。常数项记录核心交叉的带符号贡献，不是交叉根数。')+'</details>')
    out.append(heading('1. Residual from actual R3 sums','1. 实际 R3 求和得到残差',3)+'<div class="equation">R = '+esc(poly(target))+'</div>')
    out.append(para('After outside-only terms cancel, R is affine linear in signed outside-arrow counts. Its coefficients come from zero/one-outside-arrow computations. β₃ identities are interpreted in F₂.',
                    '仅选择外部箭头的项消去后，R 是外部箭头带符号计数的仿射线性函数。系数来自零／单外部箭头计算。β₃ 的恒等式在 F₂ 中解释。'))
    out.append(heading('2. Independently constructed smoothing identities','2. 独立构造平滑恒等式',3))
    out.append(para('For smoothed components A,B, H is the signed sum of crossings with A over B minus the signed sum with B over A. For classical links, both sums equal the same linking number, so H=0. This is a mathematical input, not an inference from samples.',
                    '对平滑后的分量 A、B，H 等于 A 在 B 上方时的交叉符号之和，减去 B 在 A 上方时的交叉符号之和。对经典链环，两者计算同一个 linking number，故 H=0。这是使用的数学事实，不是从样本推测出的结论。'))
    combination=[]
    for j,w in enumerate(weights):
        if not b.bv_sp.Rational(w):continue
        desc=r['identity_descriptions'][j];combination.append(f'({w}) H{j+1}')
        label=tr('({w}) H{j}: smooth crossings {smooth}; components {pair}', '({w}) H{j}：平滑交叉 {smooth}；分量 {pair}',w=esc(w),j=j+1,smooth=esc(desc['smooth']),pair=esc(desc['pair']))
        out.append('<details open><summary>'+label+'</summary>'+diagram(r['core'],components=desc['components'],smooth=desc['smooth'],spec=spec))
        out.append(para('Arc colors identify the smoothed components. Solid/dashed <b>circle arcs</b> retain the core/gap distinction. Dashed <b>chords</b> mark crossings removed by oriented smoothing. This is a Gauss-circle bookkeeping diagram, not a planar link drawing.',
                        '圆弧颜色标出平滑后的分量；实线／虚线<b>圆弧</b>仍区分核心块与补区间。虚线<b>弦</b>标出通过 oriented smoothing 去掉的交叉。这是记录分量归属的 Gauss 图，不是平面链环图。'))
        out.append('<p>'+' · '.join(f'<span style="color:{COLORS[c%len(COLORS)]}">'+tr('component {c}','分量 {c}',c=c)+'</span>' for c in sorted(set(desc['components'])))+'</p>')
        out.append('<div class="equation">H'+str(j+1)+' = '+esc(poly(r['identities'][j]))+' = 0 '+tr('<b>on classical diagrams</b>','<b>（对经典图）</b>')+'</div>')
        out.append('<p>'+' · '.join(f'B<sub>{i}</sub> → '+tr('component {c}','分量 {c}',c=c) for i,c in enumerate(desc['gap_components']))+'</p>')
        out.append('<details><summary>'+tr('Component of the segment after each endpoint','每个端点之后的线段所属分量')+'</summary><pre>'+esc(desc['components'])+'</pre></details></details>')
    out.append(heading('3. Exact coefficient comparison','3. 精确逐系数核对',3));rows=[]
    for col,expected in enumerate(target):
        actual=sum(b.bv_sp.Rational(w)*int(r['identities'][j][col]) for j,w in enumerate(weights));difference=actual-int(expected)
        assert difference%2==0 if mod else difference==0
        label=tr('constant','常数项') if col==0 else 'X['+','.join(map(str,variables[col-1]))+']'
        rows.append(f'<tr><td>{esc(label)}</td><td>{int(expected)%2 if mod else int(expected)}</td><td>{esc(actual%2 if mod else actual)}</td><td>'+tr('equal','相等')+(' (mod 2)' if mod else '')+'</td></tr>')
    out.append(para('All {n} coefficients, including the constant, are checked. This compares the complete expression; it does not require pairing individual arrows.',
                    '核对全部 {n} 个系数，包括常数项。比较的是完整表达式，不要求给每根箭头寻找配对。',n=len(target)))
    out.append('<details><summary>'+tr('Show every coefficient','展开全部系数')+'</summary><div class="table-wrap"><table><tr>'+''.join('<th>'+s+'</th>' for s in [tr('Variable','变量'),tr('R: from R3 sums','R：来自 R3 求和'),tr('Σ cᵢHᵢ: from geometry','Σ cᵢHᵢ：来自几何'),tr('Comparison','比较')])+'</tr>'+''.join(rows)+'</table></div></details>')
    out.append('<div class="equation">R = '+esc(' + '.join(combination) if combination else '0')+tr('; hence R = 0 <b>for classical completions</b>','；因此<b>对经典补全图</b>有 R = 0')+(' (mod 2)' if mod else '')+'</div>')
    return ''.join(out)


def outside_case_description(g,spec,r,name,variables,increment):
    d=outside_descriptor(g,spec);k=b.bv_names.index(name)
    column=[tuple(v) for v in variables].index((d['a'],d['b'],d['d']))+1
    coefficient=int(r['target'][column][k]);predicted=coefficient*d['sign']
    assert (int(increment)-predicted)%2==0 if name=='beta3' else int(increment)==predicted
    tail,head=(d['a'],d['b']) if d['d']==-1 else (d['b'],d['a'])
    text=tr('Outside arrow <b>{label}</b>: B<sub>{a}</sub>, B<sub>{b}</sub>; direction <b>d={d:+d}</b> (B<sub>{tail}</sub> → B<sub>{head}</sub>); crossing sign <b>{sign:+d}</b>.<br>Only X[{a},{b},{d:+d}]={sign:+d}; other outside counts are zero. Coefficient × crossing sign = ({coef}) × ({sign:+d}) = {inc}.',
            '外部箭头 <b>{label}</b>：B<sub>{a}</sub>、B<sub>{b}</sub>；方向 <b>d={d:+d}</b>（B<sub>{tail}</sub> → B<sub>{head}</sub>）；交叉符号 <b>{sign:+d}</b>。<br>仅 X[{a},{b},{d:+d}]={sign:+d}，其他外部计数为零。系数 × 交叉符号 = ({coef}) × ({sign:+d}) = {inc}。',**d,tail=tail,head=head,coef=coefficient,inc=esc(value_text(increment,name)))
    return '<p class="case-description">'+text+'</p>'


def cube_16_18_guide(r,spec):
    gg,ii=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,1],b.bv_cube.cube_table_R3_index[:,:,1],1);cards=[]
    for case,phase,term,selected,expected in [(16,1,6,[1,5],-1),(18,0,10,[3,5],1)]:
        slot=3*case+1;g=gg[:,:,phase,slot];idx=ii[phase,:,slot]
        assert int(b.bv_loop_value(gg[:,:,:,slot].transpose(2,0,1),ii[:,:,slot],[1,-1])[0])==expected
        cards.append('<div class="card"><h3>'+tr('Case {case}: β₁ residual {v:+d}','Case {case}：β₁ 残差 {v:+d}',case=case,v=expected)+'</h3>'+diagram(g,idx,selected,spec=spec)+para('p{p}, term {term}, selected arrows {arrows}.','p{p}，第 {term} 项，选中箭头 {arrows}。',p=phase+1,term=term,arrows=selected)+f'<a href="#outside-{case}">'+tr('All occurrences','查看完整贡献')+'</a></div>')
    out=heading('Why does case 16 give −1?','为什么 case 16 得到 −1？',id='case-16-18')+'<div class="grid">'+''.join(cards)+'</div>'
    out+=para('Case 16 has a negative arrow B<sub>1</sub> → B<sub>2</sub>, so X[1,2,−1]=−1. Case 18 has a negative arrow in the opposite direction, so X[1,2,+1]=−1. Their residual increments are −1 and +1.',
              'Case 16 有一根从 B<sub>1</sub> 指向 B<sub>2</sub> 的负交叉箭头，故 X[1,2,−1]=−1。Case 18 的负交叉箭头方向相反，故 X[1,2,+1]=−1。两者的残差增量分别为 −1 和 +1。')
    out+=para('The p2 value already includes path sign −1. Case 16 has only term 6: <b>(−1) × (+1) × [(+1)(−1)] × (−1) = −1</b>. Case 18 has only term 10: <b>(+1) × (+1) × [(−1)(−1)] × (+1) = +1</b>.',
              'p2 的显示值已包含路径方向 −1。Case 16 只有第 6 项：<b>(−1) × (+1) × [(+1)(−1)] × (−1) = −1</b>。Case 18 只有第 10 项：<b>(+1) × (+1) × [(−1)(−1)] × (+1) = +1</b>。')
    out+=para('If one complete diagram contains one arrow of each type, these increments cancel. <b>This does not add independent test cases as a proof, nor establish that this two-arrow completion is classical.</b>',
              '如果同一张完整图中两种箭头各有一根，它们的增量相消。<b>这不等于把独立测试 case 相加作为证明，也不证明只添加这两根箭头就能得到经典图。</b>')
    out+=para('Smooth crossings 2,3 and write A=component 1, B=component 2. In case 16, crossing 5 is the only crossing between A and B, with A over B and sign −1; in case 18, B is over A with sign −1. An isolated single crossing between two components cannot be a classical diagram. A complete classical diagram must have equal signed over/under sums.',
              '平滑交叉 2、3，记 A=分量 1，B=分量 2。Case 16 中，交叉 5 是 A、B 之间唯一的交叉，A 在 B 上方且符号为 −1；case 18 则是 B 在 A 上方且符号为 −1。这种孤立的单交叉关系不能是经典两分量图。完整经典图必须满足上下穿的带符号总和相等。')
    out+='<div class="equation">R = X[1,2,−1] − X[1,2,+1] = H7; '+tr('H7 = 0 on classical diagrams.','对经典图有 H7 = 0。')+'</div>'
    out+=para('The program compares every coefficient of R and H7; it does not find a unique partner for each arrow. The certificate below separates the geometric constraint from the exact comparison.',
              '程序逐个比较 R 与 H7 的系数，不为每根箭头寻找唯一配对。下面的证书分别展示几何约束与精确系数比较。')
    return out

@localized
def local_case(bundle,family,case,name='beta1',extra=None):
    if family=='tetra':
        if not 1<=case<=24:raise ValueError('Tetrahedron s must be 1..24.')
        r=bundle['tetra'][case-1];template=r['template'];gs=b.bv_tetra.global_tetrahedron_cases[:,:,:,template].transpose(2,0,1);ix=b.bv_tetra.global_tetrahedron_R3_index_cases[:,:,template]
        spec=arc_spec('tetra');variables=bundle['variables'];scales=[1]*8;title=f'Tetrahedron s={case} · {beta_label(name)}'
        if extra is not None:
            gg,ii=b.bv_tetra.generate_pre_tetrahedron_loops(gs.transpose(1,2,0),ix,1)
            if not 0<=extra<60:raise ValueError('Tetrahedron extra index must be 0..59.')
            gs=gg[:,:,:,extra].transpose(2,0,1);ix=ii[:,:,extra]
    elif family=='cube':
        template,rotation=case;r=next(r for r in bundle['cube'] if (r['template'],r['rotation'])==(template,rotation))
        spec=arc_spec('cube',rotation);variables=bundle['cube_variables'];scales=[1,-1];title=tr('Cube template={t}, rotation={r} · {beta}','Cube 模板={t}，基点旋转={r} · {beta}',t=template,r=rotation,beta=beta_label(name))
        if extra is None:gs=np.array([r['core'],r['core']]);ix=np.array(r['indices'])
        else:
            if not 0<=extra<40:raise ValueError('Cube outside index must be 0..39.')
            gg,ii=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,template],b.bv_cube.cube_table_R3_index[:,:,template],1)
            gs=gg[:,:,:,3*extra+rotation].transpose(2,0,1);ix=ii[:,:,3*extra+rotation]
    else:raise ValueError(family)
    if extra is not None:title+=tr(' · outside case {n}',' · 外部箭头 case {n}',n=extra)
    body=arc_guide(spec)+para('Every embedding is listed separately, including opposite contributions within a term. β₃ integer lifts are displayed; its final equation is modulo two.',
                            '每个匹配单独列出，包括同一项中的相反贡献。β₃ 展示整数提升，最终方程按模二解释。')
    if extra is not None:body+='<p class="warning">'+tr('This is a formal one-outside-arrow diagram and may be virtual. Its raw residual is retained. The certificate proves the classical equation, not vanishing on every virtual diagram.',
        '这是形式上的单外部箭头图，可能是虚拟图。其原始残差会被保留。证书证明的是经典情形的方程，不是每张虚拟图都为零。')+'</p>'
    core_residual=np.array(r['target'])[0];vals=np.array([b.bv_event(g,idx)*s for g,idx,s in zip(gs,ix,scales)]);k=b.bv_names.index(name)
    if extra is not None:body+=outside_case_description(gs[0],spec,r,name,variables,vals.sum(axis=0)[k]-core_residual[k])
    body+=heading('R3 contributions: core only' if extra is None else 'R3 contributions: selected one-arrow case','R3 贡献：无外部箭头的核心' if extra is None else 'R3 贡献：当前单箭头图')
    body+=''.join(event_report(g,idx,s,name,i+1,spec=spec) for i,(g,idx,s) in enumerate(zip(gs,ix,scales)))
    body+='<div class="equation">'+beta_label(name)+' '+tr('raw residual','原始残差')+' = '+esc(value_text(vals.sum(axis=0)[k],name))+'</div>'
    if family=='cube' and tuple(case)==(1,1) and name=='beta1' and extra is None:body+=cube_16_18_guide(r,spec)
    body+=certificate_report(r,name,variables,spec)
    if extra is None:
        body+=heading('Outside-arrow cases','逐个外部箭头 case')
        body+=para('<b>Residual increment ΔR = R(core + one outside arrow) − R(core).</b> This is not a change in the definition of β or necessarily the full residual. Each case is an independent formal input; it need not vanish individually.',
                   '<b>残差增量 ΔR = R(核心＋一根外部箭头) − R(核心)。</b>这不是 β 定义的变化，也不一定等于完整残差。每个 case 都是独立的形式输入，不要求逐个为零。')
        body+=para('{beta} core residual = {v}. Only the selected cocycle is shown.','{beta} 的核心残差 = {v}。只显示当前所选 cocycle。',beta=beta_label(name),v=esc(value_text(core_residual[k],name)))
        if family=='tetra':gg,ii=b.bv_tetra.generate_pre_tetrahedron_loops(gs.transpose(1,2,0),ix,1);slots=range(60)
        else:gg,ii=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,template],b.bv_cube.cube_table_R3_index[:,:,template],1);slots=range(rotation,120,3)
        for slot in slots:
            loop=gg[:,:,:,slot].transpose(2,0,1);indices=ii[:,:,slot];val=b.bv_loop_value(loop,indices,scales)-vals.sum(axis=0);num=slot if family=='tetra' else slot//3
            body+=f'<details id="outside-{num}"><summary>'+tr('Outside-arrow case {n} · {beta} residual increment ΔR = {v}','外部箭头 case {n} · {beta} 残差增量 ΔR = {v}',n=num,beta=beta_label(name),v=esc(value_text(val[k],name)))+'</summary>'
            body+=outside_case_description(loop[0],spec,r,name,variables,val[k])
            body+=para('Core residual {core} + increment {inc} → full one-arrow residual {total}.','核心残差 {core} ＋增量 {inc} → 完整单箭头图残差 {total}。',core=esc(value_text(core_residual[k],name)),inc=esc(value_text(val[k],name)),total=esc(value_text(core_residual[k]+val[k],name)))
            body+=''.join(event_report(g,idx,s,name,i+1,spec=spec) for i,(g,idx,s) in enumerate(zip(loop,indices,scales)))+'</details>'
    return page(title,body)

@localized
def commutation_case(u,v,a,positive=False,name=None):
    gs,ix=dr3_case(u,v,a,positive);values=np.array([b.bv_event(g,idx) for g,idx in zip(gs,ix)]);spec=arc_spec('commute')
    body=arc_guide(spec,commute=True)+para('p₁ + p₂ + p₃ + p₄ = 0. Opposite edges cancel for β₁ and β₃; β₂ may require the whole meridian.','p₁ + p₂ + p₃ + p₄ = 0。β₁、β₃ 的对边贡献相消；β₂ 可能需要完整 meridian 才相消。')
    shown=[b.bv_names.index(name)] if name else list(range(3))
    body+='<div class="table-wrap"><table><tr><th>Cocycle</th><th>p1</th><th>p2</th><th>p3</th><th>p4</th><th>'+tr('Total','总和')+'</th></tr>'
    for k in shown:
        nm=b.bv_names[k];body+='<tr><td>'+beta_label(nm)+'</td>'+''.join('<td>'+esc(value_text(v,nm))+'</td>' for v in values[:,k])+'<td>'+esc(value_text(values[:,k].sum(),nm))+'</td></tr>'
    body+='</table></div>'+''.join(event_report(g,idx,1,name,i+1,spec=spec) for i,(g,idx) in enumerate(zip(gs,ix)))
    return page(f'R3–R3 commute ({u},{v},{a})'+(tr(' · positive',' · 正交叉约定') if positive else ''),body)

@localized
def loop_case(record):
    from .loops import CFGS,TYPES,COEFF,NAMES
    body=para('The sum is calculated from the R3 events below. The predicted formula is evaluated independently using the manuscript Gauss-diagram invariants.',
              '总和来自下面各 R3 事件的实际计算。预测公式使用论文的 Gauss-diagram 不变量独立求值。')
    body+='<div class="table-wrap"><table><tr><th>Cocycle</th><th>'+tr('R3 sum','R3 总和')+'</th><th>'+tr('Predicted formula','预测公式')+'</th></tr>'+''.join('<tr><td>'+beta_label(n)+f'</td><td>{a}</td><td>'+str(p if p is not None else tr('No universal formula asserted','这里未声称一般公式'))+'</td></tr>' for n,a,p in zip(NAMES,record['value'],record['predicted']))+'</table></div>'
    body+=para('Invariants: {v}','不变量：{v}',v=esc(record['invariants']))
    if record['path_orientation']==-1:body+='<p class="warning">'+tr('The notebook rolling path is inverse to the paper convention. Raw notebook sum: {v}. Each displayed contribution includes path sign −1. β₃ is reduced modulo two.',
        '原 notebook 的 rolling 路径与论文约定相反。原始总和：{v}。显示的每个贡献已乘路径方向 −1。β₃ 最终取模二。',v=esc(record['raw_notebook_sum']))+'</p>'
    body+=para('Arrows point from underpass to overpass; arrowheads are separated from endpoint dots. Red marks the current R3 and blue the selected contributing arrows. Chord labels give crossing number and sign.',
              '箭头从下穿端指向上穿端，箭头头部与端点圆点分开。红色标出当前 R3，蓝色标出选中的贡献箭头；弦标签给出交叉编号和符号。')
    body+=heading('All R3 events','全部 R3 事件')+''.join(event_report(e['gauss'],e['indices'],record['path_orientation'],number=i+1,cfgs=CFGS,types=TYPES,coeff=COEFF,names=NAMES) for i,e in enumerate(record['events']))
    body+=heading('Complete machine-readable calculation','完整机器可读计算记录')+'<details><summary>'+tr('Gauss matrices, all occurrences, coefficients and signs','Gauss 矩阵、全部匹配、系数和符号')+'</summary>'+para('JSON field names are stable data identifiers and are shared by both language editions.','JSON 字段名是固定的数据标识符，两种语言版本共用。')+'<pre>'+esc(json.dumps(record,indent=2))+'</pre></details>'
    return page(record['kind']+' · '+record['knot']+(' · '+record['other'] if record['other'] else ''),body)

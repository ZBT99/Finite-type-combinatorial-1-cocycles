"""Actual R3 sums on push, rotation, rolling, bracket and half loops."""
import json
from pathlib import Path
import numpy as np
from . import local as b, paths, matching, invariants
from .geometry.basic_tools_jit import curl,addcurl,writhe,normalize,infty_avancer

NAMES=('alpha3','beta1','beta2','beta3')
MORTIER=np.zeros((3,4,4),dtype=np.int64)
MORTIER[:,:2,:]=np.array([[[1,1,0,0],[1,-1,0,1]],[[1,-1,0,0],[1,1,0,1]],[[1,1,0,1],[1,-1,0,3]]])
CFGS=np.concatenate((MORTIER,b.bv_cfgs));TYPES=np.r_[[4,2,2],b.bv_types]
COEFF=np.zeros((len(TYPES),4),dtype=np.int64);COEFF[:3,0]=-1;COEFF[3:,1:]=b.bv_coeff_matrix

def knots():return {n:np.array(g,dtype=np.int64) for n,g in json.loads((Path(__file__).parent/'data/knots.json').read_text()).items()}

def zero_framed(g):
    w=int(writhe(g))
    return addcurl(g,2 if w>0 else 0,abs(w)) if w else g.copy()

def expected(kind,v,u=None,w=0):
    a,c,j,e=(v[k] for k in ('v2','v3','v41','v42'));q=a*(a-1)//2
    if kind=='rotation':return [a,0,-c,0]
    if kind=='rolling':return [None,2*(q-e),2*(3*q+2*a+5*c-10*j-5*e),0]
    if kind=='half-rolling':return [None,q-e,3*q+2*a+5*c-10*j-5*e,None]
    if kind=='bracket':return [0,-2*a*u['v2'],-6*a*u['v2'],0]
    if kind=='half-bracket':return [0,-a*a,-3*a*a,0]
    if kind=='push':return [w*a,-a*u['v2'],-3*a*u['v2']-w*c,0]
    raise ValueError(kind)

def calculate(kind,name,other=None,trace=True):
    library=knots();g=library[name];h=library.get(other)
    if other=='curl+':h=curl[:,:,0].copy()
    if other=='2curls+':h=addcurl(curl[:,:,0],0,1)
    v=invariants.values(g);u=invariants.values(h) if h is not None else None
    events=[];token=matching.TRACE.set(events if trace else None);orientation=1
    try:
        if kind=='rotation':terms=paths.push_through_rl_cfgs(g,curl[:,:,0],CFGS,TYPES)
        elif kind=='rolling':
            g=zero_framed(g);terms=paths.Fox_Hatcher_rl_cfgs(g,CFGS,TYPES);orientation=-1
        elif kind=='half-rolling':
            if int(writhe(g))!=0 or not np.array_equal(normalize(g)[0],infty_avancer(g,len(g)//2)):
                raise ValueError('Half rolling requires a supplied periodic zero-framed diagram; arbitrary halving is invalid.')
            terms=paths.half_Fox_Hatcher_rl_cfgs(g,CFGS,TYPES);orientation=-1
        elif kind=='bracket':
            if h is None:raise ValueError('Bracket requires two knots.')
            g=zero_framed(g);h=zero_framed(h);terms=paths.exchange_loop_rl_cfgs(g,h,CFGS,TYPES)
        elif kind=='half-bracket':
            g=zero_framed(g);terms=paths.push_through_rl_cfgs(g,g,CFGS,TYPES)
        elif kind=='push':
            if h is None:raise ValueError('Push requires a second knot or curl.')
            terms=paths.push_through_rl_cfgs(g,h,CFGS,TYPES)
        else:raise ValueError(kind)
    finally:matching.TRACE.reset(token)
    raw=(terms@COEFF).tolist();actual=[orientation*x for x in raw];actual[3]%=2
    pred=expected(kind,v,u,int(writhe(h)) if h is not None else 0)
    if trace:
        assert np.array_equal(np.array([e['terms'] for e in events]).sum(axis=0) if events else np.zeros(len(TYPES)),terms)
        for e in events:
            e['values']=(np.array(e['terms'])@COEFF*orientation).tolist()
            for o in e['occurrences']:
                glob=o['term']-1;k=int(np.flatnonzero(COEFF[glob])[0]);o['formula']=NAMES[k]
                offset=0 if k==0 else 3+b.bv_slices[NAMES[k]].start
                o['term']=glob-offset+1;o['coefficient']=int(COEFF[glob,k])
                o['path_orientation']=orientation;o['contribution']=o['value']*o['coefficient']*orientation
    return dict(kind=kind,knot=name,other=other,invariants=v,other_invariants=u,
                gauss=g.tolist(),other_gauss=h.tolist() if h is not None else None,
                names=NAMES,raw_notebook_sum=raw,path_orientation=orientation,value=actual,predicted=pred,
                checked=[x is not None for x in pred],pass_formula=all(p is None or p==a for p,a in zip(pred,actual)),
                events=events)

def verify(output,extended=True):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    named=['3_1+','3_1-','4_1','5_1','5_2','6_1','6_2','6_3']
    examples=[]
    for n in named:
        for kind in ['rotation','rolling','half-bracket']:
            r=calculate(kind,n,trace=False);assert r['pass_formula'],r;examples.append(r)
        print('Loop checks:',n,flush=True)
    for n,m in [('4_1','3_1+'),('4_1','6_1'),('3_1+','3_1-'),('5_2','6_3')]:
        r=calculate('bracket',n,m,False);assert r['pass_formula'],r;examples.append(r)
    for n,m in [('4_1','curl+'),('4_1','2curls+'),('4_1','4_1'),('3_1+','curl+')]:
        r=calculate('push',n,m,False);assert r['pass_formula'],r;examples.append(r)
    independence=[]
    for n in ['4_1_periodic','6_1_periodic']:
        for kind in ['rotation','half-rolling']:
            r=calculate(kind,n,trace=False);assert r['pass_formula'],r;examples.append(r)
            independence.append([x%2 for x in r['value']])
    assert independence==[[1,0,0,0],[0,1,1,0],[0,0,1,0],[1,0,1,1]],independence
    # Reproduce the paper's rolling table and establish the rank used for its finite-type interpolation.
    paper_rows=[]
    for n,vals in [('3_1+',[1,1,1,0]),('4_1',[-1,0,0,0]),('5_1',[3,5,5,3]),('5_2',[2,3,3,1]),('6_1',[-2,-1,0,-1])]:
        assert list(invariants.values(knots()[n]).values())==vals
        a,c,j,e=vals;paper_rows.append([1,a,c,j,e,a*(a-1)//2])
    import sympy as sp
    assert sp.Matrix([[1,0,0,0,0,0]]+paper_rows).det()!=0
    if extended:
        for n in (f'table_{i:03}' for i in range(50)):
            for kind in ('rotation','rolling'):
                r=calculate(kind,n,trace=False);assert r['pass_formula'],r;examples.append(r)
        print('Additional knot-table checks: 100 PASS.',flush=True)
    bundle=dict(cases=len(examples),failures=0,examples=examples,independence_mod2=independence,
                interpolation_determinant=str(sp.Matrix([[1,0,0,0,0,0]]+paper_rows).det()),
                note='Loop computations test the stated examples. General formulae also use the finite-type and basis theorems of main. Raw Fox-Hatcher direction is reversed for the paper Rolling convention.')
    (output/'loop-checks.json').write_text(json.dumps(bundle,indent=2))
    return bundle

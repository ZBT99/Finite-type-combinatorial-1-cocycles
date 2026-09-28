"""The manuscript's v2, v3, v4,1=J15, v4,2=E34, by direct embeddings."""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import numpy as np
DATA=Path(__file__).parent/'data'

@lru_cache(None)
def formula(name):
    if name=='v2':return np.array([[[1,-1,0,1],[2,1,0,1],[1,1,0,1],[2,-1,0,1]]]),np.array([1])
    stem={'v3':'V3_5','v41':'J15','v42':'E34'}[name]
    ds=np.load(DATA/(stem+'.npy'));cs=np.load(DATA/(stem+'_coeff.npy'))
    if ds.shape[2]==3:ds=np.concatenate((ds,np.ones((*ds.shape[:2],1),dtype=np.int64)),axis=2)
    return ds,cs

def canonical(rows):
    mapping={};out=[]
    for a,d,*_ in rows:
        mapping.setdefault(int(a),len(mapping)+1);out.append((mapping[int(a)],int(d)))
    return tuple(out)

def evaluate(g,name,details=False):
    ds,cs=formula(name);patterns={}
    for j,(d,c) in enumerate(zip(ds,cs)):
        d=d[d[:,0]!=0];patterns.setdefault((len(d)//2,canonical(d)),[]).append((j,d,int(c)))
    labels=sorted(set(map(int,g[:,0])));out=[]
    for k in sorted({key[0] for key in patterns}):
        for selected in combinations(labels,k):
            rows=g[np.isin(g[:,0],selected)];key=(k,canonical(rows))
            for j,d,c in patterns.get(key,()):
                weight=c;seen=set();ok=True
                for (_,_,w), (a,_,req,exponent) in zip(rows,d):
                    if req and w!=req:ok=False;break
                    if a not in seen:weight*=int(w)**int(exponent);seen.add(a)
                if ok:out.append(dict(term=j+1,arrows=list(selected),coefficient=c,value=weight))
    total=sum(o['value'] for o in out)
    return (total,out) if details else total

def values(g):return {n:evaluate(g,n) for n in ('v2','v3','v41','v42')}

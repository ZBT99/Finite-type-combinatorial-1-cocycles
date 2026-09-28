"""Transparent occurrence enumeration. All arithmetic is integral."""
from itertools import combinations
from functools import lru_cache
from contextvars import ContextVar
import numpy as np
from .geometry.basic_tools_jit import get_R3_sign, r_l_R3_type

TRACE = ContextVar('trace', default=None)

def validate(g):
    g = np.asarray(g)
    if g.ndim != 2 or g.shape[1] != 3 or not np.issubdtype(g.dtype, np.integer):
        raise ValueError('A Gauss matrix must have three integer columns.')
    for label in set(g[:, 0]):
        rows=g[g[:,0]==label]
        if label<=0 or len(rows)!=2 or set(rows[:,1])!={-1,1} or rows[0,2] not in (-1,1) or rows[0,2]!=rows[1,2]:
            raise ValueError(f'Invalid crossing {label}')
    return g.astype(np.int64)

def canonical(rows, signed):
    names={}; out=[]
    for a,d,w,p in rows:
        a=int(a); names.setdefault(a,len(names)+1)
        out.append((names[a],int(d),int(w) if signed else 0,int(p)))
    return tuple(out)

@lru_cache(maxsize=32)
def _lookup(shape, data, types):
    cfgs=np.frombuffer(data,dtype=np.int64).reshape(shape)
    lookup={}
    for j,(term,typ) in enumerate(zip(cfgs,types)):
        active=term[term[:,0]!=0];signed=bool(active[0,2])
        key=(int(typ),len(active)//2,signed,canonical(active,signed))
        lookup.setdefault(key,[]).append(j)
    return lookup

def occurrences(g, index, cfgs, types):
    """Yield EVERY embedding: term index, selected labels/rows, signs, coorientation."""
    index=np.asarray(index,dtype=np.int64)
    if len(index)!=3 or not np.array_equal(index,np.sort(index)) or len(set(np.r_[index,index+1]))!=6:
        raise ValueError('R3 requires three disjoint consecutive endpoint pairs.')
    marked=set(map(int,np.r_[index,index+1])); r3=g[np.column_stack((index,index+1)).ravel()]
    typ=int(r_l_R3_type(r3));co=int(get_R3_sign(r3))
    rows=[(int(a),int(d),int(w),int(np.searchsorted(index,p,side='right'))) for p,(a,d,w) in enumerate(g) if p not in marked]
    labels=sorted({r[0] for r in rows})
    cfgs=np.asarray(cfgs,dtype=np.int64)
    lookup=_lookup(cfgs.shape,cfgs.tobytes(),tuple(map(int,types)))
    for k in sorted({key[1] for key in lookup}):
        for selected in combinations(labels,k):
            sub=[r for r in rows if r[0] in selected]
            if len(sub)!=2*k:raise ValueError('R3 does not contain complete crossings.')
            signs={r[0]:r[2] for r in sub}
            for signed in (False,True):
                weight=1 if signed else int(np.prod(list(signs.values()),dtype=np.int64))
                for j in lookup.get((typ,k,signed,canonical(sub,signed)),()):
                    yield dict(term=j+1, arrows=list(selected), rows=[list(r) for r in sub],
                               crossing_signs=list(signs.values()), signed=signed,
                               r3_type=typ, coorientation=co, matching_weight=weight,value=co*weight)

def event_terms(g,index,cfgs,types,details=False):
    occ=list(occurrences(g,index,cfgs,types));out=np.zeros(len(types),dtype=np.int64)
    for o in occ:out[o['term']-1]+=o['value']
    return (out,occ) if details else out

def R3_arr_rl_cfgs(g,index,cfgs,types):
    """Notebook-compatible R3 with a full mathematical audit trail."""
    index=np.asarray(index,dtype=np.int64)
    out=g.copy();out[index]=g[index+1];out[index+1]=g[index]
    transport=np.arange(len(g));transport[index]=index+1;transport[index+1]=index
    terms,occ=event_terms(g,index,cfgs,types,True)
    trace=TRACE.get()
    if trace is not None:
        trace.append(dict(gauss=g.tolist(),indices=index.tolist(),terms=terms.tolist(),occurrences=occ))
    return out,transport,terms

"""Reproducible finite enumeration and independently replayable exact certificates."""
import json,time
from pathlib import Path
from itertools import combinations
import numpy as np
from . import local as b

CUBE_VARIABLES=[(a,c,d) for a in range(4) for c in range(a,4) for d in (-1,1)]

def cube_certificates(details):
    records=[]
    for key,r in sorted(details.items()):
        target=np.zeros((21,3),dtype=np.int64);target[0]=r['constant']
        for j,(a,c,d) in enumerate(CUBE_VARIABLES,1):
            assert np.array_equal(r['increments'][a,c,d,-1],-r['increments'][a,c,d,1])
            target[j]=r['increments'][a,c,d,1]
        bounds=((0,2,5,8),(0,3,6,8),(0,3,5,8))[r['rotation']]
        r['target']=target
        b.bv_certify_tetra(r,CUBE_VARIABLES,[(p-1)%8 for p in bounds])
        if np.any(b.bv_in_ring(target)):
            direct=[(j,s) for j,(row,desc) in enumerate(zip(r['identities'],r['identity_descriptions']))
                    for s in (1,-1) if desc['smooth']==[2,3] and np.array_equal(s*row,target[:,0])]
            assert direct
            j,s=direct[0];weights=[0]*len(r['identities']);weights[j]=s;r['certificates']['beta1']=weights
        records.append(r)
    return records

def dr3_types(positive=False):
    table=b.bv_basic.R3_table_matrix
    types=[]
    for i in range(table.shape[2]):
        g=table[:,:,i]
        if b.bv_basic.get_R3_sign(g)!=1:continue
        if positive and not np.all(g[:,2]==1):continue
        for rotation in range(3):types.append((i,rotation,b.bv_basic.infty_avancer(g,2*rotation)))
    assert len(types)==(6 if positive else 48)
    return types

def dr3_case(u,v,arrangement,positive=False):
    types=dr3_types(positive);first,second=types[u][2],types[v][2]
    second_blocks=list(combinations(range(1,6),3))[arrangement]
    first_blocks=[a for a in range(6) if a not in second_blocks]
    p=np.array([2*a for a in first_blocks]);q=np.array([2*a for a in second_blocks])
    g=np.zeros((12,3),dtype=np.int64)
    g[np.column_stack((p,p+1)).ravel()]=first
    other=second.copy();other[:,0]+=3
    g[np.column_stack((q,q+1)).ravel()]=other
    gs=[];ix=[]
    for idx in (p,q,p,q):gs.append(g.copy());ix.append(idx);g=b.bv_deform(g,idx)
    return np.array(gs),np.array(ix)

def verify_dr3(positive=False):
    records=[];n=len(dr3_types(positive));pair_nonzero=0
    for u in range(n):
        for v in range(n):
            for a in range(10):
                gs,ix=dr3_case(u,v,a,positive);b.bv_check_closed_loop(gs,ix)
                vals=np.array([b.bv_event(g,index) for g,index in zip(gs,ix)])
                b.bv_assert_zero(vals.sum(axis=0),(u,v,a))
                pair=vals[0]+vals[2]
                assert pair[0]==0 and pair[2]%2==0,('beta1/beta3 opposite edges',(u,v,a))
                pair_nonzero+=int(pair[1]!=0)
                records.append(dict(case=[u,v,a],events=vals.tolist()))
        if u%8==0:print(f'R3-R3: first type {u+1}/{n}',flush=True)
    return dict(cases=len(records),failures=0,beta2_nonzero_opposite_edge_pairs=pair_nonzero,records=records)

def serial_record(r):
    fields=['template','rotation','core','indices','target','raw_values','identities','identity_descriptions']
    out={k:r[k].tolist() if isinstance(r[k],np.ndarray) else r[k] for k in fields if k in r}
    out['certificates']={n:[str(x) for x in r['certificates'][n]] for n in b.bv_names}
    return out

def run(output,full=True):
    start=time.time();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    cube,failures,details=b.bv_verify_cube();print('Cube enumeration:',cube,flush=True)
    cr=cube_certificates(details)
    tr=[b.bv_tetra_residual(int(t)) for t in b.bv_tetra_ids]
    for r in tr:b.bv_certify_tetra(r)
    print('Tetrahedron: all 24 cores certified.',flush=True)
    dr=verify_dr3(positive=not full)
    bundle=dict(schema='order4-certificates-v1',variables=b.bv_variables,cube_variables=CUBE_VARIABLES,
                cube_enumeration=cube,cube_raw_failures=failures,
                cube=[serial_record(r) for r in cr],tetra=[serial_record(r) for r in tr],
                dr3=dr,dr3_positive_only=not full,matcher_audit=b.bv_audit_counts,
                elapsed_seconds=round(time.time()-start,2))
    (output/'local-certificates.json').write_text(json.dumps(bundle,indent=2))
    replay(bundle)
    summary=dict(cube_cores=len(cr),cube_one_arrow_cases=5760,
                 cube_raw_nonzero=[sum(f['value'][k]%2!=0 if k==2 else f['value'][k]!=0 for f in failures) for k in range(3)],
                 tetra_cores=len(tr),tetra_one_arrow_cases=1440,
                 tetra_raw_nonzero=np.count_nonzero(b.bv_in_ring(np.concatenate([r['raw_values'] for r in tr])),axis=0).tolist(),
                 dr3_cases=dr['cases'],dr3_failures=0,certificates='PASS',matcher_audit=b.bv_audit_counts,
                 elapsed_seconds=bundle['elapsed_seconds'])
    (output/'local-summary.json').write_text(json.dumps(summary,indent=2))
    return summary

def replay(bundle):
    """Re-derive the geometric identity matrix, then compare every exact coefficient."""
    from fractions import Fraction
    for family in ('tetra','cube'):
        variables=bundle['variables' if family=='tetra' else 'cube_variables']
        for r in bundle[family]:
            core=np.array(r['core'],dtype=np.int64)
            gaps=None
            if family=='cube':
                bounds=((0,2,5,8),(0,3,6,8),(0,3,5,8))[r['rotation']];gaps=[(p-1)%8 for p in bounds]
            identities,descriptions=b.bv_linking_identities(core,variables,gaps)
            assert identities.tolist()==r['identities']
            assert descriptions==r['identity_descriptions']
            for k,name in enumerate(b.bv_names):
                w=[Fraction(x) for x in r['certificates'][name]]
                for col,expected in enumerate(np.array(r['target'])[:,k]):
                    actual=sum(w[j]*int(identities[j,col]) for j in range(len(w)))
                    delta=actual-int(expected)
                    assert delta.denominator==1 if k==2 else True
                    assert delta%2==0 if k==2 else delta==0
    return True

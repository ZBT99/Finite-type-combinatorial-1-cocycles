import unittest,json,copy
from pathlib import Path
import numpy as np
from cocycle import local as b,matching,invariants,loops
from cocycle.verification import replay,dr3_case

ROOT=Path(__file__).resolve().parents[1]

class MathematicalChecks(unittest.TestCase):
    def test_every_manuscript_vector_diagram(self):
        reference=json.loads((ROOT/'cocycle/data/manuscript-gdf.json').read_text())
        for name,prefix in [('v2','v2'),('v3','v3_'),('v41','J_'),('v42','E_')]:
            ds,cs=invariants.formula(name)
            for i,d in enumerate(ds):
                key=prefix if name=='v2' else f'{prefix}{i+1:02}'
                self.assertEqual(d[d[:,0]!=0].tolist(),reference[key]['rows'],key)
        self.assertEqual(invariants.formula('v41')[1].tolist(),[1]*10+[-1]+[1]*4)
        self.assertEqual(invariants.formula('v42')[1].tolist(),[2,2,1,-1,2,2,-1,-1,3,1,-1,1,1,3,1,1,2,-1,1,1,1,3,1,1,1,1,3,1,3,-1,1,1,-1,1])

    def test_paper_invariant_values_and_basepoint(self):
        expected={'3_1+':[1,1,1,0],'4_1':[-1,0,0,0],'5_1':[3,5,5,3],'5_2':[2,3,3,1],'6_1':[-2,-1,0,-1]}
        for name,values in expected.items():
            g=loops.knots()[name]
            for shift in range(len(g)):
                self.assertEqual(list(invariants.values(b.bv_basic.infty_avancer(g,shift)).values()),values)

    def test_detailed_embeddings_match_independent_reference(self):
        for u,v,a in [(0,0,0),(3,20,7),(47,32,9),(4,18,2)]:
            gs,ix=dr3_case(u,v,a)
            for g,index in zip(gs,ix):
                terms,occ=matching.event_terms(g,index,b.bv_cfgs,b.bv_types,True)
                np.testing.assert_array_equal(terms,b.bv_reference_event(g,index))
                np.testing.assert_array_equal(matching.event_terms(b.bv_deform(g,index),index,b.bv_cfgs,b.bv_types),-terms)
                reconstructed=np.zeros(len(terms),dtype=np.int64)
                for o in occ:reconstructed[o['term']-1]+=o['value']
                np.testing.assert_array_equal(terms,reconstructed)

    def test_virtual_cube_failure_is_not_hidden(self):
        gs,ix=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,1],b.bv_cube.cube_table_R3_index[:,:,1],1)
        bad=[]
        for s in range(gs.shape[-1]):
            delta=b.bv_event(gs[:,:,0,s],ix[0,:,s])-b.bv_event(gs[:,:,1,s],ix[1,:,s])
            if delta[0]:bad.append(delta[0])
        self.assertEqual(len(bad),4)

    def test_cube_terms_6_and_10(self):
        checked=0
        for template in (1,2,29,30):
            gs,ix=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,template],b.bv_cube.cube_table_R3_index[:,:,template],1)
            for s in range(1,gs.shape[-1],3):
                terms=[b.bv_event_terms(gs[:,:,p,s],ix[p,:,s])[b.bv_slices['beta1']]*b.bv_formulas['beta1']['coeff'] for p in (0,1)]
                delta=terms[0]-terms[1]
                self.assertEqual(int(delta.sum()),int(delta[5]+delta[9]));checked+=1
        self.assertEqual(checked,160)

    def test_replay_and_corrupted_certificate_rejection(self):
        bundle=json.loads((ROOT/'results/local-certificates.json').read_text())
        self.assertTrue(replay(bundle))
        corrupted=copy.deepcopy(bundle);corrupted['tetra'][0]['target'][0][0]+=1
        with self.assertRaises(AssertionError):replay(corrupted)
        corrupted=copy.deepcopy(bundle);corrupted['cube'][0]['identities'][0][0]+=1
        with self.assertRaises(AssertionError):replay(corrupted)

    def test_rolling_direction_and_trace_conservation(self):
        r=loops.calculate('rolling','4_1',trace=True)
        self.assertEqual(r['value'][1:],[2,2,0]);self.assertEqual(r['raw_notebook_sum'][1:3],[-2,-2])
        self.assertEqual(r['path_orientation'],-1)
        total=np.array([e['values'] for e in r['events']]).sum(axis=0);total[3]%=2
        self.assertEqual(total.tolist(),r['value'])
        contributions=[0]*4
        for e in r['events']:
            for o in e['occurrences']:contributions[loops.NAMES.index(o['formula'])]+=o['contribution']
        contributions[3]%=2;self.assertEqual(contributions,r['value'])

    def test_half_rolling_not_an_arbitrary_half(self):
        with self.assertRaises(ValueError):loops.calculate('half-rolling','3_1+')
        r=loops.calculate('half-rolling','6_1_periodic',trace=False)
        self.assertEqual([x%2 for x in r['value']],[1,0,1,1])

    def test_gauss_validation(self):
        with self.assertRaises(ValueError):matching.validate(np.array([[1,-1,1],[1,-1,1]]))
        matching.validate(loops.knots()['3_1+'])

if __name__=='__main__':unittest.main()

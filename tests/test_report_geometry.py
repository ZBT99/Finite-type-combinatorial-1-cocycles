"""Check that display labels agree with the independently stored residual coefficients."""
import json
import unittest
from pathlib import Path
import numpy as np
from cocycle import local as b
from cocycle.diagrams import arc_spec,arc_layout,outside_descriptor
from cocycle.verification import dr3_case

ROOT=Path(__file__).resolve().parents[1]

class ReportGeometryChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle=json.loads((ROOT/'results/local-certificates.json').read_text())

    def check_input(self,g,spec,record,variables,actual):
        layout=arc_layout(g,spec)
        self.assertEqual(len(layout['blocks']),len(spec['sizes']))
        self.assertEqual(len(layout['gaps']),len(spec['sizes'])+1)
        for i in range(len(g)):
            self.assertEqual(i in layout['endpoint_gaps'],int(g[i,0]) not in spec['labels'])
        desc=outside_descriptor(g,spec)
        col=[tuple(v) for v in variables].index((desc['a'],desc['b'],desc['d']))+1
        predicted=np.array(record['target'][col])*desc['sign']
        np.testing.assert_array_equal(b.bv_in_ring(actual),b.bv_in_ring(predicted))

    def test_all_cube_arc_labels_match_residual_polynomial(self):
        records={(r['template'],r['rotation']):r for r in self.bundle['cube']}
        checked=0
        for template in range(48):
            gg,ii=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,template],b.bv_cube.cube_table_R3_index[:,:,template],1)
            for slot in range(120):
                rotation=slot%3;record=records[template,rotation];spec=arc_spec('cube',rotation)
                loop=gg[:,:,:,slot].transpose(2,0,1)
                actual=b.bv_loop_value(loop,ii[:,:,slot],[1,-1])-np.array(record['target'][0])
                self.check_input(loop[0],spec,record,self.bundle['cube_variables'],actual);checked+=1
        self.assertEqual(checked,5760)

    def test_all_tetra_arc_labels_match_residual_polynomial(self):
        checked=0
        for record in self.bundle['tetra']:
            t=record['template'];base=b.bv_tetra.global_tetrahedron_cases[:,:,:,t];ix=b.bv_tetra.global_tetrahedron_R3_index_cases[:,:,t]
            gg,ii=b.bv_tetra.generate_pre_tetrahedron_loops(base,ix,1)
            for slot in range(60):
                loop=gg[:,:,:,slot].transpose(2,0,1)
                for g in loop:arc_layout(g,arc_spec('tetra'))
                actual=b.bv_loop_value(loop,ii[:,:,slot])-np.array(record['target'][0])
                self.check_input(loop[0],arc_spec('tetra'),record,self.bundle['variables'],actual);checked+=1
        self.assertEqual(checked,1440)

    def test_commute_six_blocks_survive_each_event(self):
        for u,v in ((0,2),(47,32)):
            for arrangement in range(10):
                gs,ix=dr3_case(u,v,arrangement)
                for g,index in zip(gs,ix):
                    layout=arc_layout(g,arc_spec('commute'))
                    self.assertEqual(layout['blocks'],[(0,1),(2,3),(4,5),(6,7),(8,9),(10,11)])
                    for start in index:self.assertIn((int(start),int(start)+1),layout['blocks'])

    def test_case_16_and_18_smoothing_explanation(self):
        gg,ii=b.bv_cube.generate_pre_cube_loops(b.bv_cube.cube_table_matrix[:,:,1],b.bv_cube.cube_table_R3_index[:,:,1],1)
        for case,direction,expected in [(16,-1,-1),(18,1,1)]:
            slot=3*case+1;g=gg[:,:,0,slot]
            desc=outside_descriptor(g,arc_spec('cube',1))
            self.assertEqual(desc,dict(label=5,a=1,b=2,d=direction,sign=-1))
            component,count=b.bv_oriented_smoothing(g,[2,3]);self.assertEqual(count,3)
            between=[]
            for label in (1,4,5):
                ends=np.flatnonzero(g[:,0]==label)
                if set(map(int,component[ends]))=={1,2}:between.append(label)
            self.assertEqual(between,[5])
            self.assertEqual(int(b.bv_loop_value(gg[:,:,:,slot].transpose(2,0,1),ii[:,:,slot],[1,-1])[0]),expected)

if __name__=='__main__':unittest.main()

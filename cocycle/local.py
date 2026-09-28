if not __debug__:
    raise RuntimeError("Exact verification must not be run with Python -O (assertions disabled).")

import sys as bv_sys
import json as bv_json
import hashlib as bv_hashlib
import itertools as bv_it
from pathlib import Path as bv_Path
import numpy as bv_np
import sympy as bv_sp
bv_display = print
bv_Markdown = str

bv_repo = bv_Path(__file__).resolve().parent
bv_formula_dir = bv_repo / 'data'
from .geometry import basic_tools_jit as bv_basic
from .geometry import gen_cube_loop as bv_cube
from .geometry import gen_tetra_loop as bv_tetra
try:
    import os
    if os.environ.get('COCYCLE_NATIVE_AUDIT') == '0':
        raise ImportError('Native audit disabled for portability check')
    import basic_tools_c as bv_btc
except ImportError:
    bv_btc = None  # The independent Python matcher below also runs on its own.
bv_code_hashes = {bv_Path(m.__file__).name: bv_hashlib.sha256(bv_Path(m.__file__).read_bytes()).hexdigest()
                  for m in (bv_basic, bv_cube, bv_tetra) if m.__file__}

bv_formulas = {}
bv_source_hashes = {}
for bv_name, bv_stem, bv_subdir, bv_modulus in (
    ('beta1', 'beta1', '', None),
    ('beta2', 'beta2', '', None),
    ('beta3', 'Z2', 'modp', 2),
):
    bv_arrays = []
    for bv_prefix in ('arr_cfgs', 'arr_R3_types', 'arr_coeff'):
        bv_path = bv_formula_dir / bv_subdir / f'{bv_prefix}_{bv_stem}.npy'
        bv_array = bv_np.load(bv_path, allow_pickle=False)
        assert bv_np.issubdtype(bv_array.dtype, bv_np.integer), bv_path
        bv_arrays.append(bv_np.ascontiguousarray(bv_array, dtype=bv_np.int64))
        bv_source_hashes[bv_path.name] = bv_hashlib.sha256(bv_path.read_bytes()).hexdigest()
    bv_cfg, bv_types, bv_coeff = bv_arrays
    assert bv_cfg.shape == (len(bv_coeff), 4, 4)
    assert bv_types.shape == bv_coeff.shape
    assert bv_np.all((0 <= bv_types) & (bv_types < 6))
    for bv_term in bv_cfg:
        bv_active = bv_term[bv_term[:, 0] != 0]
        assert len(bv_active) in (2, 4)
        assert (len(bv_active) == 4 and bv_np.all(bv_active[:, 2] == 0)) or (
            len(bv_active) == 2 and bv_np.all(bv_np.abs(bv_active[:, 2]) == 1)
        ), 'The finite outside-arrow reduction below assumes unsigned2 + signed1.'
    bv_formulas[bv_name] = dict(cfg=bv_cfg, types=bv_types, coeff=bv_coeff, modulus=bv_modulus)

bv_names = tuple(bv_formulas)
bv_cfgs = bv_np.concatenate([bv_formulas[n]['cfg'] for n in bv_names])
bv_types = bv_np.concatenate([bv_formulas[n]['types'] for n in bv_names])
bv_coeff_matrix = bv_np.zeros((len(bv_types), 3), dtype=bv_np.int64)
bv_slices = {}
bv_offset = 0
for bv_k, bv_name in enumerate(bv_names):
    bv_end = bv_offset + len(bv_formulas[bv_name]['coeff'])
    bv_slices[bv_name] = slice(bv_offset, bv_end)
    bv_coeff_matrix[bv_offset:bv_end, bv_k] = bv_formulas[bv_name]['coeff']
    bv_offset = bv_end

def bv_canonical_rows(rows, signed):
    labels, result = {}, []
    for label, direction, sign, piece in rows:
        label = int(label)
        labels.setdefault(label, len(labels) + 1)
        result.append((labels[label], int(direction), int(sign) if signed else 0, int(piece)))
    return tuple(result)

bv_lookup = {}
for bv_j, (bv_term, bv_type) in enumerate(zip(bv_cfgs, bv_types)):
    bv_active = bv_term[bv_term[:, 0] != 0]
    bv_signed = bool(bv_active[0, 2] != 0)
    bv_key = (int(bv_type), len(bv_active)//2, bv_signed, bv_canonical_rows(bv_active, bv_signed))
    bv_lookup.setdefault(bv_key, []).append(bv_j)

def bv_reference_event(g, index):
    index = bv_np.asarray(index, dtype=bv_np.int64)
    assert bv_np.array_equal(index, bv_np.sort(index))
    marked = bv_np.column_stack((index, index+1)).ravel()
    assert len(set(map(int, marked))) == 6
    r3 = g[marked]
    assert len(set(map(int, r3[:, 0]))) == 3
    typ = int(bv_basic.r_l_R3_type(r3))
    coorientation = int(bv_basic.get_R3_sign(r3))
    rows = [(int(a), int(d), int(w), int(bv_np.searchsorted(index, p, side='right')))
            for p, (a, d, w) in enumerate(g) if p not in marked]
    labels = sorted({r[0] for r in rows})
    result = bv_np.zeros(len(bv_cfgs), dtype=bv_np.int64)
    for k in (1, 2):
        for selected in bv_it.combinations(labels, k):
            sub = [r for r in rows if r[0] in selected]
            assert len(sub) == 2*k
            signs = {r[0]: r[2] for r in sub}
            unsigned_weight = coorientation
            for w in signs.values():
                unsigned_weight *= w
            for signed in (False, True):
                key = (typ, k, signed, bv_canonical_rows(sub, signed))
                for j in bv_lookup.get(key, ()):
                    result[j] += coorientation if signed else unsigned_weight
    return result

bv_audit_counts = {'events': 0, 'native_comparisons': 0}

def bv_event_terms(g, index):
    g = bv_np.ascontiguousarray(g, dtype=bv_np.int64)
    index = bv_np.ascontiguousarray(index, dtype=bv_np.int64)
    reference = bv_reference_event(g, index)
    bv_audit_counts['events'] += 1
    if bv_btc is not None:
        native = bv_btc.R3_arr_rl_cfgs(g, index, bv_cfgs, bv_types)[2]
        assert bv_np.array_equal(native, reference), 'Python/C term mismatch'
        bv_audit_counts['native_comparisons'] += 1
    return reference

def bv_event(g, index):
    return bv_event_terms(g, index) @ bv_coeff_matrix

def bv_deform(g, index):
    out = g.copy()
    out[index] = g[index+1]
    out[index+1] = g[index]
    return out

def bv_check_closed_loop(cases, indices):
    for phase in range(len(cases)):
        assert bv_np.array_equal(
            bv_basic.normalize(bv_deform(cases[phase], indices[phase]))[0],
            bv_basic.normalize(cases[(phase+1) % len(cases)])[0]
        ), ('The event sequence is not a closed R3 loop', phase)

def bv_loop_value(cases, indices, scales=None):
    if scales is None:
        scales = [1]*len(cases)
    return sum((int(s)*bv_event(g, idx) for g, idx, s in zip(cases, indices, scales)),
               start=bv_np.zeros(3, dtype=bv_np.int64))

def bv_in_ring(values):
    out = bv_np.array(values, dtype=bv_np.int64, copy=True)
    out[..., 2] %= 2
    return out

def bv_assert_zero(values, context):
    reduced = bv_in_ring(values)
    assert not bv_np.any(reduced), (context, 'integer residual', values.tolist(),
                                   'residual in the correct rings', reduced.tolist())

def bv_verify_cube():
    records, failures, details = [], [], {}
    for extra in (0, 1):
        count = 0
        for template in range(bv_cube.cube_table_matrix.shape[2]):
            gs, ix = bv_cube.generate_pre_cube_loops(
                bv_cube.cube_table_matrix[:, :, template],
                bv_cube.cube_table_R3_index[:, :, template], extra
            )
            for case in range(gs.shape[-1]):
                g = gs[:, :, 0, case]
                # Both branches start at the same diagram.
                assert bv_np.array_equal(g, gs[:, :, 1, case])
                indices = ix[:, :, case]
                signatures = []
                for idx in indices:
                    r3 = g[bv_np.column_stack((idx, idx+1)).ravel()]
                    outside_pieces = tuple(int(bv_np.searchsorted(idx, p, side='right'))
                                           for p in bv_np.flatnonzero(g[:, 0] > 4))
                    signatures.append((int(bv_basic.r_l_R3_type(r3)),
                                       int(bv_basic.get_R3_sign(r3)), outside_pieces))
                assert signatures[0] == signatures[1]
                value = bv_event(g, indices[0]) - bv_event(g, indices[1])
                rotation = case % 3
                key = (template, rotation)
                detail = details.setdefault(key, dict(template=template, rotation=rotation, increments={}))
                if extra == 0:
                    detail.update(core=g.copy(), indices=indices.copy(), constant=value.copy())
                else:
                    p, q = bv_np.flatnonzero(g[:, 0] == 5)
                    boundaries = ((0, 2, 5, 8), (0, 3, 6, 8), (0, 3, 5, 8))[rotation]
                    a, b = boundaries.index(int(p)), boundaries.index(int(q-1))
                    direction, sign = map(int, g[p, 1:3])
                    detail['increments'][a, b, direction, sign] = value-detail['constant']
                if bv_np.any(bv_in_ring(value)):
                    failures.append(dict(template=template, case=case, extra=extra, value=value.tolist()))
                count += 1
        records.append(dict(equation='cube', extra=extra, cases=count, failures=len(failures)))
    assert [r['cases'] for r in records] == [144, 5760]
    return records, failures, details

def bv_positive_dr3_cases():
    local_positive = [i for i in range(bv_basic.R3_table_matrix.shape[2])
                      if bv_np.all(bv_basic.R3_table_matrix[:, 2, i] == 1)
                      and bv_basic.get_R3_sign(bv_basic.R3_table_matrix[:, :, i]) == 1]
    types = [bv_basic.infty_avancer(bv_basic.R3_table_matrix[:, :, i], 2*r)
             for i in local_positive for r in range(3)]
    assert len(types) == 6
    assert {int(bv_basic.r_l_R3_type(g)) for g in types} == set(range(6))
    for u, first in enumerate(types):
        for v, second in enumerate(types):
            for arrangement, second_blocks in enumerate(bv_it.combinations(range(1, 6), 3)):
                first_blocks = [a for a in range(6) if a not in second_blocks]
                p = bv_np.array([2*a for a in first_blocks], dtype=bv_np.int64)
                q = bv_np.array([2*a for a in second_blocks], dtype=bv_np.int64)
                g = bv_np.zeros((12, 3), dtype=bv_np.int64)
                g[bv_np.column_stack((p, p+1)).ravel()] = first
                other = second.copy()
                other[:, 0] += 3
                g[bv_np.column_stack((q, q+1)).ravel()] = other
                cases, indices = [], []
                for idx in (p, q, p, q):
                    cases.append(g.copy())
                    indices.append(idx.copy())
                    g = bv_deform(g, idx)
                yield (u, v, arrangement), bv_np.array(cases), bv_np.array(indices)

def bv_verify_dr3():
    failures, active = [], [0, 0, 0]
    count = 0
    for label, gs, ix in bv_positive_dr3_cases():
        bv_check_closed_loop(gs, ix)
        values = bv_np.array([bv_event(g, idx) for g, idx in zip(gs, ix)])
        for k in range(3):
            active[k] += int(bv_np.any(values[:, k] % 2 if k == 2 else values[:, k]))
        if bv_np.any(bv_in_ring(values.sum(axis=0))):
            failures.append(dict(case=label, events=values.tolist()))
        count += 1
    assert count == 360
    assert not failures, failures[:5]
    return dict(equation='positive DR3', extra=0, cases=count, failures=0,
                cases_with_nonzero_event_values=dict(zip(bv_names, active)))

bv_tetra_ids = bv_np.flatnonzero(bv_np.all(bv_tetra.global_tetrahedron_cases[:, 2, :, :] == 1,
                                        axis=(0, 1)))
assert len(bv_tetra_ids) == 24
bv_variables = [(a, b, d) for a in range(5) for b in range(a, 5) for d in (-1, 1)]
bv_variable_index = {key: j for j, key in enumerate(bv_variables)}

def bv_tetra_outside_cancellation(loop, indices):
    signatures = []
    for g, idx in zip(loop, indices):
        r3 = g[bv_np.column_stack((idx, idx+1)).ravel()]
        # Insert a hypothetical endpoint immediately before each local 3-endpoint block,
        # or after the final block. No crossing endpoint lies inside an R3 pair.
        pieces = tuple(int(bv_np.count_nonzero(idx < gap)) for gap in (0, 3, 6, 9, 12))
        signatures.append((int(bv_basic.r_l_R3_type(r3)), int(bv_basic.get_R3_sign(r3)), pieces))
    for j in range(4):
        a, b = signatures[j], signatures[j+4]
        assert a[0] == b[0] and a[1] == -b[1] and a[2] == b[2]

def bv_tetra_residual(template):
    g = bv_tetra.global_tetrahedron_cases[:, :, :, template]
    ix = bv_tetra.global_tetrahedron_R3_index_cases[:, :, template]
    loop, indices = g.transpose(2, 0, 1), ix
    bv_check_closed_loop(loop, indices)
    bv_tetra_outside_cancellation(loop, indices)
    constant = bv_loop_value(loop, indices)
    gs, extra_ix = bv_tetra.generate_pre_tetrahedron_loops(g, ix, 1)
    assert gs.shape[-1] == 60
    increments, raw_values = {}, []
    for case in range(gs.shape[-1]):
        extra_loop = gs[:, :, :, case].transpose(2, 0, 1)
        bv_check_closed_loop(extra_loop, extra_ix[:, :, case])
        p, q = bv_np.flatnonzero(extra_loop[0, :, 0] == 7)
        key = (int(p//3), int((q-1)//3), int(extra_loop[0, p, 1]), int(extra_loop[0, p, 2]))
        assert key not in increments
        residual = bv_loop_value(extra_loop, extra_ix[:, :, case])
        increments[key] = residual - constant
        raw_values.append(residual)
    target = bv_np.zeros((31, 3), dtype=bv_np.int64)
    target[0] = constant
    for j, (a, b, d) in enumerate(bv_variables, 1):
        assert bv_np.array_equal(increments[a, b, d, -1], -increments[a, b, d, 1])
        target[j] = increments[a, b, d, 1]
    return dict(template=int(template), core=loop[0].copy(), loop=loop.copy(),
                indices=indices.copy(), target=target, raw_values=bv_np.array(raw_values))


def bv_oriented_smoothing(core, smooth):
    length = len(core)
    successor = (bv_np.arange(length) + 1) % length
    for label in smooth:
        p, q = bv_np.flatnonzero(core[:, 0] == label)
        successor[(p-1) % length] = q
        successor[(q-1) % length] = p
    assert sorted(map(int, successor)) == list(range(length))
    component = bv_np.full(length, -1, dtype=bv_np.int64)
    count = 0
    for segment in range(length):
        if component[segment] >= 0:
            continue
        p = segment
        while component[p] < 0:
            component[p] = count
            p = successor[p]
        assert p == segment
        count += 1
    return component, count

def bv_linking_identities(core, variables=None, gap_segments=None):
    if variables is None:
        variables = bv_variables
    if gap_segments is None:
        gap_segments = [11, 2, 5, 8, 11]
    identities = {}
    labels = sorted(set(map(int, core[:, 0])))
    # Prefer fewer smoothings, so the resulting certificates are easy to inspect.
    for number in range(1, len(labels)+1):
        for smooth in bv_it.combinations(labels, number):
            component, count = bv_oriented_smoothing(core, smooth)
            for a, b in bv_it.combinations(range(count), 2):
                def crossing_factor(p, q, direction, sign):
                    if {int(component[p]), int(component[q])} != {a, b}:
                        return 0
                    # +1 is the arrow head (overpassing endpoint) in these Gauss matrices.
                    return int((1 if component[p] == a else -1)*direction*sign)
                constant = 0
                for label in labels:
                    if label not in smooth:
                        p, q = bv_np.flatnonzero(core[:, 0] == label)
                        constant += crossing_factor(p, q, core[p, 1], core[p, 2])
                row = [constant] + [crossing_factor(gap_segments[l], gap_segments[r], d, 1)
                                    for l, r, d in variables]
                if not any(row):
                    continue
                key = tuple(row)
                identities.setdefault(key, dict(smooth=list(map(int, smooth)), pair=[a, b],
                                                components=component.tolist(),
                                                gap_components=component[gap_segments].tolist()))
    assert identities
    return bv_np.array(list(identities), dtype=bv_np.int64), list(identities.values())

def bv_solve_gf2(matrix, rhs):
    a = bv_np.column_stack((matrix, rhs)).astype(bv_np.int64) % 2
    row, pivots = 0, []
    for col in range(matrix.shape[1]):
        nonzero = bv_np.flatnonzero(a[row:, col])
        if not len(nonzero):
            continue
        pivot = row + nonzero[0]
        a[[row, pivot]] = a[[pivot, row]]
        for j in range(len(a)):
            if j != row and a[j, col]:
                a[j] ^= a[row]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    assert not bv_np.any(a[row:, -1]), 'No linking identity certificate over F_2'
    solution = bv_np.zeros(matrix.shape[1], dtype=bv_np.int64)
    solution[pivots] = a[:row, -1]
    assert not bv_np.any((matrix @ solution-rhs) % 2)
    return solution.tolist()

def bv_solve_exact(matrix, rhs):
    matrix, rhs = bv_sp.Matrix(matrix), bv_sp.Matrix(rhs)
    solution, parameters = matrix.gauss_jordan_solve(rhs)
    solution = solution.subs({p: 0 for p in parameters})
    assert matrix*solution == rhs
    return list(solution)

def bv_certify_tetra(record, variables=None, gap_segments=None):
    h, descriptions = bv_linking_identities(record['core'], variables, gap_segments)
    certificates = {}
    for k, name in enumerate(bv_names):
        target = record['target'][:, k]
        weights = bv_solve_gf2(h.T, target) if k == 2 else bv_solve_exact(h.T, target)
        certificates[name] = weights
        # Reconstruct every coefficient, including the constant, without any tolerance.
        reconstructed = [sum(weights[j]*int(h[j, col]) for j in range(len(h)))
                         for col in range(h.shape[1])]
        differences = [value-int(expected) for value, expected in zip(reconstructed, target)]
        assert all(v % 2 == 0 if k == 2 else v == 0 for v in differences)
    record.update(identities=h, identity_descriptions=descriptions, certificates=certificates)
    return record



def bv_polynomial_text(vector, modulus=None, variables=None):
    if variables is None:
        variables = bv_variables
    parts = []
    for j, value in enumerate(vector):
        value = bv_sp.Rational(value)
        if modulus:
            value = int(value) % modulus
        if not value:
            continue
        if j == 0:
            term = '1'
        else:
            a, b, d = variables[j-1]
            term = f'X[{a},{b},{d:+d}]'
        parts.append(f'({value})*{term}')
    return ' + '.join(parts) if parts else '0'

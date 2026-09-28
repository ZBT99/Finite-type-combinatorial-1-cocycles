"""Display-only Gauss layouts; block boundaries come from the local generators."""
import html
import itertools
import math
import numpy as np
from .language import tr

CUBE_BLOCKS = ((2, 3, 3), (3, 3, 2), (3, 2, 3))
COLORS = ('#176c9c', '#b94c12', '#248044', '#8050a0', '#ab465b', '#5c8080')
_DIAGRAM_IDS = itertools.count()


def arc_spec(family, rotation=0):
    if family == 'cube':
        return dict(sizes=CUBE_BLOCKS[rotation], labels=tuple(range(1, 5)))
    if family == 'tetra':
        return dict(sizes=(3, 3, 3, 3), labels=tuple(range(1, 7)))
    if family == 'commute':
        return dict(sizes=(2,) * 6, labels=tuple(range(1, 7)))
    raise ValueError(family)


def arc_layout(g, spec):
    """Classify based gaps without confusing them with a single R3's four intervals."""
    g = np.asarray(g)
    core = [i for i, row in enumerate(g) if int(row[0]) in spec['labels']]
    if len(core) != sum(spec['sizes']):
        raise ValueError('Wrong number of core endpoints for the requested blocks')
    blocks = []
    offset = 0
    for size in spec['sizes']:
        positions = core[offset:offset + size]
        if positions != list(range(positions[0], positions[-1] + 1)):
            raise ValueError('An outside endpoint lies inside a fixed core block')
        blocks.append((positions[0], positions[-1]))
        offset += size
    gaps = [(-1, blocks[0][0])]
    gaps += [(left[1], right[0]) for left, right in zip(blocks, blocks[1:])]
    gaps += [(blocks[-1][1], len(g))]
    endpoint_gaps = {}
    for b, (left, right) in enumerate(gaps):
        for i in range(left + 1, right):
            endpoint_gaps[i] = b
    # Fixed gap widths keep B_0 and the last B separate and legible around infinity.
    weights = [1.] * (len(g) + 1)
    segment_gaps = {}
    for b, (left, right) in enumerate(gaps):
        width = 1.6 if b in (0, len(gaps) - 1) else 2.3
        for segment in range(left + 1, right + 1):
            weights[segment] = width / (right - left)
            segment_gaps[segment] = b
    total = sum(weights)
    boundaries = [0.]
    for weight in weights:
        boundaries.append(boundaries[-1] + 2 * math.pi * weight / total)
    return dict(blocks=blocks, gaps=gaps, endpoint_gaps=endpoint_gaps,
                angles=boundaries[1:-1], boundaries=boundaries, segment_gaps=segment_gaps)


def outside_descriptor(g, spec):
    """The actual signed unit input to one X[a,b,d], extracted from a case matrix."""
    g = np.asarray(g)
    layout = arc_layout(g, spec)
    labels = sorted(set(map(int, g[:, 0])) - set(spec['labels']))
    if len(labels) != 1:
        raise ValueError('Expected exactly one outside arrow')
    label = labels[0]
    p, q = map(int, np.flatnonzero(g[:, 0] == label))
    return dict(label=label, a=layout['endpoint_gaps'][p], b=layout['endpoint_gaps'][q],
                d=int(g[p, 1]), sign=int(g[p, 2]))


def diagram(g, index=(), selected=None, components=None, smooth=(), spec=None):
    """Solid core blocks, dashed based gaps; chord arrows always point under to over."""
    g = np.asarray(g)
    n = len(g)
    cx = cy = 210
    radius = 146
    layout = arc_layout(g, spec) if spec else None
    angles = layout['angles'] if layout else [2 * math.pi * (i + .5) / n for i in range(n)]

    def point(angle, r=radius):
        return cx - r * math.sin(angle), cy - r * math.cos(angle)

    def arc(start, end, color, dashed=False, kind='', number=None):
        x, y = point(start)
        u, v = point(end)
        data = f' data-arc="{kind}"'
        if number is not None:
            data += f' data-gap="{number}"'
        dash = ' stroke-dasharray="5 5"' if dashed else ''
        return (f'<path{data} d="M{x:.2f},{y:.2f} A{radius},{radius} 0 '
                f'{int(end-start > math.pi)} 0 {u:.2f},{v:.2f}" fill="none" '
                f'stroke="{color}" stroke-width="2.5"{dash}/>')

    marked = {int(g[i, 0]) for p in index for i in (p, p + 1)}
    selected = set(selected or ())
    marker_id = f'arrowhead-{next(_DIAGRAM_IDS)}'
    pts = [point(a) for a in angles]
    description = tr('Core blocks are solid arcs; B gaps are dashed. Core endpoints are filled; outside endpoints are hollow. Read counterclockwise from infinity. Arrows point from underpass to overpass.', '核心块为实线弧，B 补区间为虚线弧。核心端点实心，外部端点空心。从基点逆时针读取；箭头从下穿端指向上穿端。')
    out = ['<svg class="gauss-diagram" viewBox="0 0 420 420" role="img" '
           f'aria-label="{html.escape(description)}"><title>{tr("Based Gauss diagram", "带基点的 Gauss 图")}</title>'
           f'<defs><marker id="{marker_id}" viewBox="0 0 10 10" refX="9" refY="5" '
           'markerUnits="userSpaceOnUse" markerWidth="11" markerHeight="9" orient="auto-start-reverse">'
           '<path d="M0,0 L10,5 L0,10z" fill="context-stroke"/></marker></defs>']
    if layout:
        boundaries = layout['boundaries']
        for segment in range(n + 1):
            gap = layout['segment_gaps'].get(segment)
            color = COLORS[int(components[(segment - 1) % n]) % len(COLORS)] if components is not None else '#465466'
            out.append(arc(boundaries[segment], boundaries[segment + 1], color,
                           gap is not None, 'gap' if gap is not None else 'core', gap))
        for b, (left, right) in enumerate(layout['gaps']):
            start = 0 if left < 0 else angles[left]
            end = 2 * math.pi if right == n else angles[right]
            x, y = point((start + end) / 2, radius + 38)
            out.append(f'<text class="gap-label" x="{x:.2f}" y="{y+5:.2f}" '
                       f'text-anchor="middle" fill="#176c9c" font-size="16" font-weight="bold">'
                       f'B<tspan baseline-shift="sub" font-size="12">{b}</tspan></text>')
    else:
        out.append(f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="#465466"/>')
    for label in sorted(set(g[:, 0])) if n else ():
        ends = np.flatnonzero(g[:, 0] == label)
        p = int(next(i for i in ends if g[i, 1] == -1))
        q = int(next(i for i in ends if g[i, 1] == 1))
        x, y = pts[p]
        u, v = pts[q]
        color = '#ba3333' if label in marked else '#167cad' if label in selected else '#53606e'
        if selected and label not in selected and label not in marked:
            color = '#7c8794'
        # Stop the shaft before the endpoint dot, so its arrowhead stays visible.
        length = math.hypot(u-x, v-y)
        ux, uy = (u-x)/length, (v-y)/length
        tip_x, tip_y = u-7*ux, v-7*uy
        dashed = ' stroke-dasharray="3 3"' if label in smooth else ''
        out.append(f'<path data-arrow="{int(label)}" d="M{x:.2f},{y:.2f} L{tip_x:.2f},{tip_y:.2f}" '
                   f'fill="none" stroke="{color}" stroke-width="{2.3 if label in selected or label in marked else 1.5}" '
                   f'marker-end="url(#{marker_id})"{dashed}/>')
        # Put arrow numbers close to their tails, away from central chord intersections.
        out.append(f'<text x="{.82*x+.18*u+4:.2f}" y="{.82*y+.18*v-5:.2f}" '
                   f'fill="{color}" font-size="14">{int(label)}:{"+" if g[p,2]>0 else "−" if g[p,2]<0 else "±"}</text>')
    for i, ((x, y), angle) in enumerate(zip(pts, angles)):
        is_core = not spec or int(g[i, 0]) in spec['labels']
        out.append(f'<circle data-endpoint="{i}" cx="{x:.2f}" cy="{y:.2f}" r="3.4" '
                   f'fill="{"#465466" if is_core else "white"}" stroke="#465466" stroke-width="1.2"/>')
        u, v = point(angle, radius + 13)
        out.append(f'<text x="{u:.2f}" y="{v+4:.2f}" text-anchor="middle" '
                   f'fill="#667" font-size="11">{i}</text>')
    out.append(f'<circle cx="{cx}" cy="{cy-radius}" r="4" fill="#176c9c"/>'
               f'<text x="{cx}" y="{cy-radius-17}" text-anchor="middle" font-size="18" fill="#176c9c">∞</text></svg>')
    return ''.join(out)

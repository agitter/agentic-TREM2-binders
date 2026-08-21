#!/usr/bin/env python3
"""
Combined network: all 100 TREM2 designs in one layout.

Same layout strategy and same edge encoding as fig_network.py; the only change is
that node hue now carries the wet-lab outcome (binder vs non-binder) instead of the
design workflow, so within- and across-category fold similarity is visible in one
picture. Unlike the split panels, this figure plots the binder<->non-binder edges.

Usage:  python3 fig_network_combined.py [TM_THR] [ID_THR]
"""
import sys
import numpy as np, pandas as pd, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

import fig_network as F
from render_nodes import render, superpose, CA, _spline
from py2Dmol.viewer import best_view

DPI = F.DPI
plt.rcParams.update({'savefig.bbox': 'tight', 'lines.scale_dashes': False})

# node hue = outcome. Chosen for maximum hue separation from each other (146 deg)
# and from the three edge colours; see METHODS_ARCHIVE.md 7c.
NODE = {True: '#4A3F9E', False: '#E08A1E'}
NODE_LABEL = {True: 'binder', False: 'non-binder'}

REF = F.order[int(np.argmax(F.TM.sum(1)))]
FRAME = best_view(_spline(CA[REF]))

_cache = {}
def node_image(d):
    if d not in _cache:
        C = superpose(REF, d)
        if C is None:
            C = CA[d]
        _cache[d] = render(C, NODE[F.status[d]], px=300, width=12.0, frame=FRAME)
    return _cache[d]


def pair_class(u, v):
    a, b = F.status[u], F.status[v]
    return 'bb' if (a and b) else ('nn' if not (a or b) else 'bn')


def draw(nodes, pos, EG, title, subtitle, statline, statline2, fname,
         min_d, width_in=17.0, accent='#2B3A4A',
         wmin=0.0, wmax=1.0, lw_base=2.1, lw_gain=6.6):
    xs = np.array([pos[n][0] for n in nodes]); ys = np.array([pos[n][1] for n in nodes])
    pad = min_d * 0.62
    spanx = (xs.max() - xs.min()) + 2 * pad
    spany = (ys.max() - ys.min()) + 2 * pad
    AXF = 0.720
    fig = plt.figure(figsize=(width_in, width_in * spany / spanx / AXF))
    ax = fig.add_axes([0.0, 0.150, 1.0, AXF])

    for u, v, dd in sorted(EG.edges(data=True), key=lambda e: e[2]['w']):
        a = np.clip((dd['w'] - wmin) / (wmax - wmin + 1e-9), 0, 1)
        k = F.edge_class(u, v)
        af = {'diff': 0.68, 'same': 0.90, 'unknown': 0.90}[k]
        ax.plot(*zip(pos[u], pos[v]), color=F.EDGE_COLORS[k],
                lw=lw_base + lw_gain * a ** 1.3,
                alpha=af + (0.98 - af) * a ** 1.1,
                zorder=(1 if k == 'diff' else 2) + a, solid_capstyle='round')

    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    ax.set_ylim(ys.min() - pad, ys.max() + pad)
    ax.set_aspect('equal'); ax.set_axis_off()

    fig.canvas.draw()
    p0 = ax.transData.transform((0, 0)); p1 = ax.transData.transform((min_d / 2, 0))
    img_zoom = (2 * abs(p1[0] - p0[0])) / 300 * (72.0 / DPI) * 1.16
    for n in nodes:
        ax.add_artist(AnnotationBbox(OffsetImage(node_image(n), zoom=img_zoom),
                                     pos[n], frameon=False, zorder=3))

    fig.suptitle(title, fontsize=40, weight='bold', y=0.998)
    fig.text(.5, .955, subtitle, ha='center', fontsize=22, color='#3A3A3A')
    fig.text(.5, .928, statline, ha='center', fontsize=24, color=accent, weight='bold')
    if statline2:
        fig.text(.5, .897, statline2, ha='center', fontsize=19, color='#4A4A4A')

    n_b = sum(1 for n in nodes if F.status[n])
    h1 = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=NODE[True],
                     markersize=22, label=f'binder  ({n_b})'),
          plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=NODE[False],
                     markersize=22, label=f'non-binder  ({len(nodes) - n_b})')]
    leg1 = ax.legend(handles=h1, loc='lower center', bbox_to_anchor=(.5, -.055), ncol=2,
                     frameon=False, fontsize=22, handletextpad=.3, columnspacing=2.4)
    ax.add_artist(leg1)

    counts = {k: 0 for k in F.EDGE_COLORS}
    for u, v in EG.edges():
        counts[F.edge_class(u, v)] += 1
    h2 = [plt.Line2D([0], [0], color=F.EDGE_COLORS[k], lw=6.0,
                     label=f'{F.EDGE_LABEL[k]}  ({counts[k]})')
          for k in ('diff', 'same', 'unknown') if counts[k]]
    ax.legend(handles=h2, loc='lower center', bbox_to_anchor=(.5, -.122),
              ncol=len(h2), frameon=False, fontsize=19, handlelength=3.0,
              handletextpad=.5, columnspacing=2.2)

    plt.savefig(fname, dpi=DPI, facecolor='white')
    plt.close()
    print(f'{fname}: {len(nodes)} nodes, {EG.number_of_edges()} edges')


if __name__ == '__main__':
    TM_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.80
    ID_THR = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
    TAG = f'_tm{int(round(TM_THR * 100)):02d}'

    nodes = list(F.order)                       # all 100, binders and non-binders
    SG = F.build(nodes, TM_THR)
    QG = F.build(nodes, ID_THR, M=F.SID)
    MIN_D = 0.098                               # 100 nodes rather than 81
    pos, cid, comms = F.force_layout(SG, nodes, MIN_D, seed=7)

    def split(G):
        c = {'bb': 0, 'bn': 0, 'nn': 0}
        for u, v in G.edges():
            c[pair_class(u, v)] += 1
        return c

    cs, cq = split(SG), split(QG)
    both = sum(1 for u, v in SG.edges() if QG.has_edge(u, v))
    indep = [(u, v) for u, v in SG.edges() if F.edge_class(u, v) != 'same']
    indep_both = sum(1 for u, v in indep if QG.has_edge(u, v))

    draw(nodes, pos, SG,
         'Fold space of all TREM2 designs',
         'each node is a designed protein, drawn from its predicted structure  ·  '
         'all 100 designs, binders and non-binders together',
         f'{SG.number_of_edges()} pairs share a fold  (TM-score ≥ {TM_THR:.2f})',
         f'{cs["bb"]} binder–binder  ·  {cs["bn"]} binder–non-binder  ·  '
         f'{cs["nn"]} non-binder–non-binder',
         f'fig_network_combined{TAG}_structure.png',
         min_d=MIN_D, accent='#2B3A4A', wmin=TM_THR, wmax=1.0)

    draw(nodes, pos, QG,
         'Sequence space of all TREM2 designs',
         'same designs, same positions  ·  lines now join pairs whose '
         'amino-acid sequences are alike',
         f'only {both} of those {SG.number_of_edges()} fold-sharing pairs are also '
         f'alike in sequence  (identity ≥ {ID_THR:.0f}%)',
         f'{cq["bb"]} binder–binder  ·  {cq["bn"]} binder–non-binder  ·  '
         f'{cq["nn"]} non-binder–non-binder  ·  '
         f'only {indep_both} of {len(indep)} independent pairs overlap',
         f'fig_network_combined{TAG}_sequence.png',
         min_d=MIN_D, accent='#B23A22', wmin=ID_THR, wmax=100.0)

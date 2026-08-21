#!/usr/bin/env python3
"""
Three-collection fold-space and sequence-space networks (200 designs).

Node hue = collection.  Edge colour = backbone provenance, exactly as in
fig_network.py; the hackathon and Muni collections publish no backbone
identifiers, so any edge touching them is 'backbone not recorded'.

Usage: python3 fig_network_hack.py [TM_THR] [ID_THR]
"""
import sys, json, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from py2Dmol.viewer import best_view

import fig_network as F                      # layout engine
import render_nodes as R

# Same three edge classes as fig_network.py, but the balance is inverted here:
# only Anthropic publishes backbone ids, so 'not recorded' is now the majority
# class rather than a rare special case.  It is therefore drawn as a light
# recessive neutral so the two informative classes stay legible.
EDGE_COLORS = {'diff': '#2A3340', 'same': '#D81B60', 'unknown': '#8FD4D8'}
EDGE_ALPHA  = {'diff': 0.80, 'same': 0.92, 'unknown': 0.42}

DPI = 125
plt.rcParams.update({'savefig.bbox': 'tight', 'lines.scale_dashes': False})

TM = np.load('tm_hack.npy'); SID = np.load('sid_hack.npy')
order = json.load(open('order_hack.json'))
idx = {n: i for i, n in enumerate(order)}
prov = pd.read_csv('provenance_summary.csv')
BACKBONE = dict(zip(prov.full_name, prov.root_backbone_id))

NODE = {'Hackathon': '#1F9E8C', 'Muni': '#C1440E', 'Anthropic': '#4A3F9E'}
LABEL = {'Hackathon': 'Adaptyv × MUNI hackathon', 'Muni': 'Muni autoresearch',
         'Anthropic': 'Anthropic (Claude)'}


def collection(d):
    if d.startswith('HACK_'):
        return 'Hackathon'
    if d.startswith('MUNI_'):
        return 'Muni'
    return 'Anthropic'


def edge_class(u, v):
    """Only Anthropic designs carry backbone ids, so every other edge is unknown."""
    if collection(u) != 'Anthropic' or collection(v) != 'Anthropic':
        return 'unknown'
    return 'same' if BACKBONE[u] == BACKBONE[v] else 'diff'


def build(nodes, thr, M):
    G = nx.Graph(); G.add_nodes_from(nodes)
    for i, u in enumerate(nodes):
        for v in nodes[i + 1:]:
            t = float(M[idx[u], idx[v]])
            if t >= thr:
                G.add_edge(u, v, w=t)
    return G


# ---- node bitmaps, shared orientation frame, coloured by collection
REF = order[int(np.argmax(TM.sum(1)))]
_ca_cache = {}


def _load_ca():
    """Cα traces + Foldseek CIGARs for the 200-design set."""
    import os, glob
    CA = {}
    for f in glob.glob('all_binder_pdb_hack/*.pdb'):
        xs = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
              for l in open(f) if l.startswith('ATOM') and l[12:16].strip() == 'CA']
        CA[os.path.basename(f)[:-4]] = np.array(xs)
    aln = pd.read_csv('fs_hack_sup/sup.tsv', sep='\t', header=None,
                      names=['q', 't', 'tm', 'qs', 'qe', 'ts', 'te', 'cigar'])
    aln['key'] = aln.q + '|' + aln.t
    return CA, aln.set_index('key')


CA, ALN = _load_ca()
FRAME = best_view(R._spline(CA[REF]))


def superpose(ref, mob):
    if ref == mob:
        return CA[ref].copy()
    k = ref + '|' + mob
    if k not in ALN.index:
        return None
    r = ALN.loc[k]
    r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
    qi, ti, P = int(r.qs) - 1, int(r.ts) - 1, []
    for n, op in re.findall(r'(\d+)([MID])', r.cigar):
        n = int(n)
        if op == 'M':
            P += [(qi + k2, ti + k2) for k2 in range(n)]; qi += n; ti += n
        elif op == 'I':
            qi += n
        else:
            ti += n
    if len(P) < 4:
        return None
    X = CA[ref][[i for i, _ in P]]; Y = CA[mob][[j for _, j in P]]
    H = (Y - Y.mean(0)).T @ (X - X.mean(0))
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (CA[mob] - Y.mean(0)) @ Rm.T + X.mean(0)


def node_image(d):
    if d not in _ca_cache:
        C = superpose(REF, d)
        if C is None:
            C = CA[d]
        _ca_cache[d] = R.render(C, NODE[collection(d)], px=300, width=12.0, frame=FRAME)
    return _ca_cache[d]


def draw(nodes, pos, EG, title, subtitle, statline, statline2, fname,
         min_d, width_in=17.0, accent='#2B3A4A', wmin=0.0, wmax=1.0,
         lw_base=2.1, lw_gain=6.6):
    xs = np.array([pos[n][0] for n in nodes]); ys = np.array([pos[n][1] for n in nodes])
    pad = min_d * 0.62
    spanx = (xs.max() - xs.min()) + 2 * pad
    spany = (ys.max() - ys.min()) + 2 * pad
    AXF = 0.720
    fig = plt.figure(figsize=(width_in, width_in * spany / spanx / AXF))
    ax = fig.add_axes([0.0, 0.150, 1.0, AXF])

    for u, v, dd in sorted(EG.edges(data=True), key=lambda e: e[2]['w']):
        a = np.clip((dd['w'] - wmin) / (wmax - wmin + 1e-9), 0, 1)
        k = edge_class(u, v)
        af = EDGE_ALPHA[k]
        ax.plot(*zip(pos[u], pos[v]), color=EDGE_COLORS[k],
                lw=(lw_base + lw_gain * a ** 1.3) * (0.72 if k == 'unknown' else 1.0),
                alpha=af + (0.98 - af) * a ** 1.1,
                zorder=(0 if k == 'unknown' else 2) + a, solid_capstyle='round')

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

    cnt = {c: sum(1 for n in nodes if collection(n) == c) for c in NODE}
    h1 = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=NODE[c],
                     markersize=22, label=f'{LABEL[c]}  ({cnt[c]})') for c in NODE]
    leg1 = ax.legend(handles=h1, loc='lower center', bbox_to_anchor=(.5, -.055), ncol=3,
                     frameon=False, fontsize=20, handletextpad=.3, columnspacing=1.6)
    ax.add_artist(leg1)
    counts = {k: 0 for k in F.EDGE_COLORS}
    for u, v in EG.edges():
        counts[edge_class(u, v)] += 1
    h2 = [plt.Line2D([0], [0], color=EDGE_COLORS[k], lw=6.0,
                     label=f'{F.EDGE_LABEL[k]}  ({counts[k]})')
          for k in ('diff', 'same', 'unknown') if counts[k]]
    ax.legend(handles=h2, loc='lower center', bbox_to_anchor=(.5, -.122),
              ncol=len(h2), frameon=False, fontsize=17, handlelength=3.0,
              handletextpad=.45, columnspacing=1.8)
    plt.savefig(fname, dpi=DPI, facecolor='white'); plt.close()
    print(f'{fname}: {len(nodes)} nodes, {EG.number_of_edges()} edges')


if __name__ == '__main__':
    TM_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.80
    ID_THR = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
    TAG = f'_tm{int(round(TM_THR * 100)):02d}'
    MIN_D = 0.072                                  # 200 nodes

    SG = build(order, TM_THR, TM)
    QG = build(order, ID_THR, SID)
    pos, cid, comms = F.force_layout(SG, order, MIN_D, seed=7)

    def cross(G):
        c = collections_counter(G)
        return c
    def collections_counter(G):
        c = {}
        for u, v in G.edges():
            a, b = sorted([collection(u), collection(v)])
            c[f'{a}-{b}'] = c.get(f'{a}-{b}', 0) + 1
        return c

    cs, cq = collections_counter(SG), collections_counter(QG)
    def fmt(c, keys):
        return '  ·  '.join(f'{k.replace("Hackathon","Hack")}: {c.get(k,0)}' for k in keys)
    keys = ['Anthropic-Hackathon', 'Hackathon-Muni', 'Anthropic-Muni']

    draw(order, pos, SG,
         'Fold space of 200 TREM2 designs',
         'three independent campaigns, one layout  ·  '
         'nodes drawn from their Boltz-2 co-folded structures',
         f'{SG.number_of_edges()} pairs share a fold  (TM-score ≥ {TM_THR:.2f})',
         'cross-campaign links — ' + fmt(cs, keys),
         f'fig_network_hack{TAG}_structure.png',
         min_d=MIN_D, width_in=20.0, accent='#2B3A4A', wmin=TM_THR, wmax=1.0)

    both = sum(1 for u, v in SG.edges() if QG.has_edge(u, v))
    draw(order, pos, QG,
         'Sequence space of the same 200 designs',
         'same designs, same positions  ·  lines now join pairs whose '
         'amino-acid sequences are alike',
         f'only {both} of those {SG.number_of_edges()} fold-sharing pairs are also '
         f'alike in sequence  (identity ≥ {ID_THR:.0f}%)',
         'cross-campaign links — ' + fmt(cq, keys),
         f'fig_network_hack{TAG}_sequence.png',
         min_d=MIN_D, width_in=20.0, accent='#B23A22', wmin=ID_THR, wmax=100.0)

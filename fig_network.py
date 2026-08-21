import numpy as np, pandas as pd, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from render_nodes import render, superpose, CA, _spline
from py2Dmol.viewer import best_view

DPI = 125
# Matplotlib scales dash patterns by linewidth, so a thick dotted line renders with
# gaps several times longer and reads as dashed. Disable that: line style must encode
# backbone provenance only, independent of the linewidth that encodes similarity.
plt.rcParams.update({'savefig.bbox': 'tight', 'lines.scale_dashes': False})

TM = np.load('tm_all.npy')
SID = np.load('sid_all.npy')          # pairwise % sequence identity, same order
order = json.load(open('order_all.json'))
idx = {n: i for i, n in enumerate(order)}
arch = pd.read_csv('archetypes.csv').set_index('design')
_prov = pd.read_csv('provenance_summary.csv')
BACKBONE = dict(zip(_prov.full_name, _prov.root_backbone_id))


def edge_class(u, v):
    """Is this edge independent evidence, or two sequences on one backbone?

    Muni carries no backbone identifiers, so a Muni-Muni pair is genuinely unknown.
    A Muni-Anthropic pair cannot share a backbone (separate pipelines), so it counts
    as different-backbone."""
    mu_u, mu_v = u.startswith('MUNI_'), v.startswith('MUNI_')
    if mu_u and mu_v:
        return 'unknown'
    if mu_u != mu_v:
        return 'diff'
    return 'same' if BACKBONE[u] == BACKBONE[v] else 'diff'


# Edge COLOUR encodes backbone provenance; every edge is solid so connectivity stays
# legible. Chosen to avoid the four node hues. See METHODS_ARCHIVE.md 7c.
EDGE_COLORS = {'diff': '#2A3340', 'same': '#D81B60', 'unknown': '#00C2C7'}
EDGE_LABEL = {'diff': 'different backbones', 'same': 'same backbone',
              'unknown': 'backbone not recorded'}
mu = pd.read_json('muni.json').set_index('id')

PAL = {'muni': '#C1440E', 'mythos_preview_single': '#0E7C7B',
       'mythos_preview_multi': '#2E5EAA', 'opus_4_8_multi': '#8E5AA8'}
LB = {'muni': 'Muni', 'mythos_preview_single': 'Mythos-Preview single-target',
      'mythos_preview_multi': 'Mythos-Preview multi-target', 'opus_4_8_multi': 'Opus-4.8 multi-target'}

status = {}
for d in order:
    if d.startswith('MUNI_'):
        status[d] = (mu.loc[d[5:], 'binding_strength'] == 'Strong')
    else:
        status[d] = bool(arch.loc[d, 'binder'])

REF = order[int(np.argmax(TM.sum(1)))]
FRAME = best_view(_spline(CA[REF]))

_cache = {}
def node_image(d):
    if d not in _cache:
        C = superpose(REF, d)
        if C is None:
            C = CA[d]
        _cache[d] = render(C, PAL[arch.loc[d, 'workflow']], px=300, width=12.0, frame=FRAME)
    return _cache[d]


def build(nodes, thr, M=None):
    M = TM if M is None else M
    G = nx.Graph(); G.add_nodes_from(nodes)
    for i, u in enumerate(nodes):
        for v in nodes[i + 1:]:
            t = float(M[idx[u], idx[v]])
            if t >= thr:
                G.add_edge(u, v, w=t)
    return G


def force_layout(G, nodes, min_d, seed=7, iters=1500):
    """Community-seeded force layout.

    Edges attract, everything repels, and overlaps are resolved by the same
    physics rather than a separate pass, so module structure survives.
    """
    rng = np.random.default_rng(seed)
    comms = nx.community.louvain_communities(G, weight='w', seed=seed, resolution=1.0)
    # sets iterate in hash order, which varies per process; impose a total order
    comms = [sorted(c) for c in comms]
    comms.sort(key=lambda c: (-len(c), c[0]))
    cid = {n: k for k, c in enumerate(comms) for n in c}

    # lay out the community graph first, so modules get their own territory
    CG = nx.Graph(); CG.add_nodes_from(range(len(comms)))
    for u, v, d in G.edges(data=True):
        if cid[u] != cid[v]:
            a, b = cid[u], cid[v]
            CG.add_edge(a, b, w=CG.get_edge_data(a, b, {'w': 0})['w'] + d['w'])
    cpos = nx.spring_layout(CG, weight='w', seed=seed, k=1.6, iterations=600) \
        if len(comms) > 1 else {0: np.zeros(2)}

    P = np.zeros((len(nodes), 2))
    for i, n in enumerate(nodes):
        P[i] = np.asarray(cpos[cid[n]], float) * 1.5 + rng.normal(0, 0.09, 2)

    ni = {n: i for i, n in enumerate(nodes)}
    E = np.array([[ni[u], ni[v]] for u, v in G.edges()], dtype=int).reshape(-1, 2)
    W = np.array([G[u][v]['w'] for u, v in G.edges()], dtype=float)
    if len(E):
        W = (W - W.min()) / (np.ptp(W) + 1e-9) * 0.8 + 0.4
    L = min_d * 1.02

    for it in range(iters):
        cool = 1.0 - 0.75 * it / iters
        F = np.zeros_like(P)
        if len(E):
            d = P[E[:, 1]] - P[E[:, 0]]
            dist = np.hypot(d[:, 0], d[:, 1])[:, None] + 1e-9
            f = (dist - L) / L * W[:, None] * 0.045
            np.add.at(F, E[:, 0], f * d / dist)
            np.add.at(F, E[:, 1], -f * d / dist)
        diff = P[:, None, :] - P[None, :, :]
        dist = np.hypot(diff[..., 0], diff[..., 1])
        np.fill_diagonal(dist, 9e9)
        u = diff / (dist[..., None] + 1e-9)
        rep = np.clip((min_d * 1.7) / (dist + 1e-9) - 1.0, 0, None)[..., None] * 0.0045
        F += (u * rep).sum(1)
        F += (u * np.clip(min_d - dist, 0, None)[..., None] * 0.55).sum(1)
        F -= P * 0.0030
        P += np.clip(F, -min_d * 0.5, min_d * 0.5) * cool

    # pack disconnected components: lay each out on its own, then place the
    # blobs by decreasing size so the page stays compact instead of sprawling
    comps = [sorted(c) for c in nx.connected_components(G)]
    comps.sort(key=lambda c: (-len(c), c[0]))
    if len(comps) > 1:
        blobs = []
        for comp in comps:
            ii = [ni[n] for n in comp]
            Q = P[ii] - P[ii].mean(0)
            r = np.hypot(Q[:, 0], Q[:, 1]).max() + min_d * 0.55
            blobs.append((list(comp), Q, r))
        placed = []
        for k, (members, Q, r) in enumerate(blobs):
            if k == 0:
                c = np.zeros(2)
            else:
                best, bestcost = None, 9e18
                for ang in np.linspace(0, 2 * np.pi, 180, endpoint=False):
                    for rad in np.arange(0.05, 14, 0.05):
                        c = np.array([np.cos(ang), np.sin(ang)]) * rad
                        ok = all(np.hypot(*(c - pc)) >= (r + pr) * 1.02
                                 for pc, pr in placed)
                        if ok:
                            cost = np.hypot(*c)
                            if cost < bestcost:
                                bestcost, best = cost, c
                            break
                c = best if best is not None else np.array([0., 0.])
            placed.append((c, r))
            for n, q in zip(members, Q):
                P[ni[n]] = c + q

    # settle: gravity pulls the satellites inward until they touch the core,
    # edge springs hold each module's internal shape, collisions keep it legible
    for it in range(900):
        F = np.zeros_like(P)
        if len(E):
            d = P[E[:, 1]] - P[E[:, 0]]
            dist = np.hypot(d[:, 0], d[:, 1])[:, None] + 1e-9
            f = (dist - L) / L * W[:, None] * 0.055
            np.add.at(F, E[:, 0], f * d / dist)
            np.add.at(F, E[:, 1], -f * d / dist)
        diff = P[:, None, :] - P[None, :, :]
        dist = np.hypot(diff[..., 0], diff[..., 1]); np.fill_diagonal(dist, 9e9)
        u = diff / (dist[..., None] + 1e-9)
        F += (u * np.clip(min_d - dist, 0, None)[..., None] * 0.60).sum(1)
        F -= (P - P.mean(0)) * 0.010
        P += np.clip(F, -min_d * 0.35, min_d * 0.35)

    for _ in range(800):
        diff = P[:, None, :] - P[None, :, :]
        dist = np.hypot(diff[..., 0], diff[..., 1]); np.fill_diagonal(dist, 9e9)
        if not (dist < min_d).any():
            break
        u = diff / (dist[..., None] + 1e-9)
        P += (u * np.clip(min_d - dist, 0, None)[..., None] * 0.5).sum(1)
    return {n: P[i] for i, n in enumerate(nodes)}, cid, comms


def draw(nodes, pos, EG, title, subtitle, statline, fname,
         min_d=0.118, width_in=17.0, edge_color='#39404B',
         wmin=0.0, wmax=1.0, lw_base=2.1, lw_gain=6.6, alpha_base=0.58,
         statline2=None):
    """Draw one panel. Node positions are passed in, so the structural and the
    sequence panel of a set are the same picture with a different edge layer."""
    xs = np.array([pos[n][0] for n in nodes]); ys = np.array([pos[n][1] for n in nodes])
    pad = min_d * 0.62
    spanx = (xs.max() - xs.min()) + 2 * pad
    spany = (ys.max() - ys.min()) + 2 * pad
    AXF = 0.735
    fig = plt.figure(figsize=(width_in, width_in * spany / spanx / AXF))
    ax = fig.add_axes([0.0, 0.150, 1.0, AXF])

    for u, v, dd in sorted(EG.edges(data=True), key=lambda e: e[2]['w']):
        a = np.clip((dd['w'] - wmin) / (wmax - wmin + 1e-9), 0, 1)
        k = edge_class(u, v)
        af = {'diff': 0.68, 'same': 0.90, 'unknown': 0.90}[k]
        ax.plot(*zip(pos[u], pos[v]), color=EDGE_COLORS[k],
                lw=lw_base + lw_gain * a ** 1.3,
                alpha=af + (0.98 - af) * a ** 1.1,
                zorder=(1 if k == 'diff' else 2) + a,
                solid_capstyle='round')

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
    if statline:
        fig.text(.5, .928, statline, ha='center', fontsize=24,
                 color=edge_color, weight='bold')
    if statline2:
        fig.text(.5, .897, statline2, ha='center', fontsize=19, color='#4A4A4A')
    h = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                    markersize=22, label=LB[k]) for k, c in PAL.items()]
    leg1 = ax.legend(handles=h, loc='lower center', bbox_to_anchor=(.5, -.055), ncol=4,
                     frameon=False, fontsize=20, handletextpad=.3, columnspacing=1.3)
    ax.add_artist(leg1)
    counts = {k: 0 for k in EDGE_COLORS}
    for u, v in EG.edges():
        counts[edge_class(u, v)] += 1
    h2 = [plt.Line2D([0], [0], color=EDGE_COLORS[k], lw=6.0,
                     label=f'{EDGE_LABEL[k]}  ({counts[k]})')
          for k in ('diff', 'same', 'unknown') if counts[k]]
    ax.legend(handles=h2, loc='lower center', bbox_to_anchor=(.5, -.122),
              ncol=len(h2), frameon=False, fontsize=19, handlelength=3.0,
              handletextpad=.5, columnspacing=2.2)
    plt.savefig(fname, dpi=DPI, facecolor='white')
    plt.close()
    print(f'{fname}: {len(nodes)} nodes, {EG.number_of_edges()} edges')


if __name__ == '__main__':
    import sys
    TM_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.80
    ID_THR = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
    TAG = f'_tm{int(round(TM_THR * 100)):02d}'

    SETS = [
        dict(nodes=[d for d in order if status[d]], key='binders', min_d=0.108,
             what='binders',
             noun='the 81 designs that bind TREM2'),
        dict(nodes=[d for d in order if not status[d]], key='nonbinders', min_d=0.215,
             what='non-binders',
             noun='the 19 designs that failed to bind TREM2'),
    ]

    for S in SETS:
        nodes, min_d = S['nodes'], S['min_d']
        SG = build(nodes, TM_THR)                       # structural graph
        QG = build(nodes, ID_THR, M=SID)                # sequence graph
        pos, cid, comms = force_layout(SG, nodes, min_d, seed=7)

        n_s, n_q = SG.number_of_edges(), QG.number_of_edges()
        both = sum(1 for u, v in SG.edges() if QG.has_edge(u, v))
        sib_s = sum(1 for u, v in SG.edges() if edge_class(u, v) == 'same')
        indep_s = n_s - sib_s
        indep_both = sum(1 for u, v in SG.edges()
                         if QG.has_edge(u, v) and edge_class(u, v) != 'same')

        draw(nodes, pos, SG,
             f'Fold space of TREM2 {S["what"]}',
             f'each node is a designed protein, drawn from its predicted structure  ·  '
             f'{S["noun"]}',
             f'{n_s} pairs share a fold  (TM-score ≥ {TM_THR:.2f})',
             f'fig_network_{S["key"]}{TAG}_structure.png',
             min_d=min_d, edge_color='#2B3A4A', wmin=TM_THR, wmax=1.0,
             statline2=f'{indep_s} join independent designs; {sib_s} join two sequences '
                       f'built on one backbone')

        draw(nodes, pos, QG,
             f'Sequence space of the same TREM2 {S["what"]}',
             f'same designs, same positions  ·  lines now join pairs whose '
             f'amino-acid sequences are alike',
             f'only {both} of those {n_s} fold-sharing pairs are also alike in sequence  '
             f'(identity ≥ {ID_THR:.0f}%)',
             f'fig_network_{S["key"]}{TAG}_sequence.png',
             min_d=min_d, edge_color='#B23A22', wmin=ID_THR, wmax=100.0,
             lw_base=2.1, lw_gain=6.6, alpha_base=0.58,
             statline2=f'and {both - indep_both} of those {both} are sibling sequences on a '
                       f'shared backbone — only {indep_both} of {indep_s} independent pairs')

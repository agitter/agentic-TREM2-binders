#!/usr/bin/env python3
"""Structure (TM) and sequence (% identity) similarity heatmaps, 200 designs.

Two orderings are produced:
  _bycollection : blocked by campaign, hierarchically clustered within each block
  _clustered    : one global UPGMA ordering, campaign shown only as a side bar
"""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

plt.rcParams.update({'savefig.bbox': 'tight'})
COL = {'Hackathon': '#1F9E8C', 'Muni': '#C1440E', 'Anthropic': '#4A3F9E'}
LBL = {'Hackathon': 'Adaptyv × MUNI hackathon', 'Muni': 'Muni autoresearch',
       'Anthropic': 'Anthropic (Claude)'}

TM = np.load('tm_hack.npy'); SID = np.load('sid_hack.npy')
order = json.load(open('order_hack.json')); idx = {n: i for i, n in enumerate(order)}
h = pd.read_csv('hackathon_designs.csv')
CLS = dict(zip('HACK_' + h.pb_id, h.design_class))


def coll(d):
    return 'Hackathon' if d.startswith('HACK_') else ('Muni' if d.startswith('MUNI_') else 'Anthropic')


def is_ab(d):
    return CLS.get(d) in ('Nanobody', 'scFv')


def upgma(M, ii):
    if len(ii) < 3:
        return list(ii)
    D = 1 - M[np.ix_(ii, ii)] if M.max() <= 1.001 else (100 - M[np.ix_(ii, ii)]) / 100
    np.fill_diagonal(D, 0); D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), 'average')
    return [ii[k] for k in dendrogram(Z, no_plot=True)['leaves']]


def draw(fname, ordering, blocks, title_note):
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 8.6))
    for ax, M, cmap, vmin, vmax, lab, ttl in [
            (axes[0], TM, 'viridis', 0.2, 1.0, 'TM-score',
             'Structural similarity  (Foldseek TM-align)'),
            (axes[1], SID, 'magma', 0, 60, '% sequence identity',
             'Sequence similarity  (BLOSUM62, shorter-sequence normalised)')]:
        im = ax.imshow(M[np.ix_(ordering, ordering)], cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(ttl, fontsize=13.5, weight='bold', pad=8)
        n = len(ordering)
        for k, i in enumerate(ordering):                       # campaign side bar
            ax.add_patch(plt.Rectangle((-6.5, k - .5), 5, 1, color=COL[coll(order[i])],
                                       clip_on=False))
        for (s, e, name) in blocks:                            # campaign separators
            ax.axhline(e - .5, color='white', lw=1.4)
            ax.axvline(e - .5, color='white', lw=1.4)
        # the 5 antibody-scaffold designs are outlined only; no text label, so the
        # matrix itself is never covered.  Explained in the sub-title instead.
        abpos = [k for k, i in enumerate(ordering) if is_ab(order[i])]
        if abpos and max(abpos) - min(abpos) == len(abpos) - 1:
            a, b = min(abpos) - .5, max(abpos) + .5
            ax.add_patch(plt.Rectangle((a, a), b - a, b - a, fill=False,
                                       ec='#FF2D95', lw=2.4, zorder=5))
        cb = plt.colorbar(im, ax=ax, fraction=.046, pad=.02)
        cb.set_label(lab, fontsize=12); cb.ax.tick_params(labelsize=10)
    hs = [plt.Line2D([0], [0], color=COL[c], lw=7, label=f'{LBL[c]}') for c in COL]
    fig.legend(handles=hs, loc='lower center', ncol=3, frameon=False, fontsize=12.5,
               bbox_to_anchor=(.5, -.035))
    fig.suptitle('TREM2 designs: structure vs sequence similarity across three campaigns',
                 fontsize=19, weight='bold', y=1.015)
    fig.text(.5, .965, title_note, ha='center', fontsize=12.5, color='#444')
    plt.tight_layout()
    plt.savefig(fname, dpi=150, facecolor='white'); plt.close()
    print(fname)


if __name__ == '__main__':
    # ---- ordering 1: blocked by campaign, clustered within
    ordering, blocks, cur = [], [], 0
    for c in ['Hackathon', 'Muni', 'Anthropic']:
        ii = [i for i, n in enumerate(order) if coll(n) == c]
        # antibodies first inside the hackathon block so they form one contiguous run
        if c == 'Hackathon':
            ab = [i for i in ii if is_ab(order[i])]
            rest = [i for i in ii if not is_ab(order[i])]
            ii = upgma(TM, ab) + upgma(TM, rest)
        else:
            ii = upgma(TM, ii)
        ordering += ii; blocks.append((cur, cur + len(ii), c)); cur += len(ii)
    draw('fig_heatmap_hack_bycollection.png', ordering, blocks,
         'ordered by campaign, clustered within each block  ·  pink outline = the 5 '
         'antibody-scaffold designs (max TM 0.58 to any other design)')

    # ---- ordering 2: one global clustering
    g = upgma(TM, list(range(len(order))))
    draw('fig_heatmap_hack_clustered.png', g, [],
         'single global UPGMA ordering on TM-score, campaign shown only in the side bar  ·  '
         'pink outline = the 5 antibody-scaffold designs')

#!/usr/bin/env python3
"""Novelty audit (Anthropic and Muni, separately) + three-collection epitope map."""
import json, collections, itertools
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'savefig.bbox': 'tight'})
COL = {'Hackathon': '#1F9E8C', 'Muni': '#C1440E', 'Anthropic': '#4A3F9E'}

TM = np.load('tm_hack.npy'); SID = np.load('sid_hack.npy')
order = json.load(open('order_hack.json')); idx = {n: i for i, n in enumerate(order)}


# ------------------------------------------------------------------ figure 1
def novelty_audit():
    d = pd.read_csv('novelty_audit.csv')
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.0))
    panels = [('Anthropic', 'Anthropic designs (n = 90)',
               'the novelty gate in Anthropic’s prompt'),
              ('Muni', 'Muni autoresearch designs (n = 10)',
               'Muni published no equivalent gate — shown for comparison')]
    for ax, (tag, title, note) in zip(axes, panels):
        g = d[d.collection == tag]
        ax.axhspan(0.80, 1.02, color='#F3D9DE', zorder=0)
        ax.axvspan(60, 100, color='#E8E8E8', zorder=0)
        ax.axhline(0.50, ls='--', lw=1.2, color='#555')
        ax.axhline(0.80, ls='-', lw=1.4, color='#C0392B')
        ax.axvline(60, ls='-', lw=1.8, color='#333')
        ax.scatter(g.best_sid_global, g.best_tm, s=70, c=COL[tag],
                   edgecolor='white', lw=.7, zorder=3)
        ax.set_xlim(0, 100); ax.set_ylim(0.30, 1.02)
        ax.set_xlabel('best sequence identity to any hackathon design (%)', fontsize=13)
        ax.set_ylabel('best TM-score to any hackathon design', fontsize=13)
        n80 = int((g.best_tm >= 0.80).sum())
        ax.set_title(f'{title}\n{n80} of {len(g)} are structurally near-identical '
                     f'to a public design\nyet none is rejected by sequence novelty',
                     fontsize=13, weight='bold')
        ax.text(59, 0.325, 'novelty gate rejects →\n(>60% id, >50% cov)', fontsize=10,
                ha='right', va='bottom', color='#333')
        ax.text(1.5, 0.81, 'TM ≥ 0.80  near-identical backbone', fontsize=10,
                color='#C0392B', va='bottom')
        ax.text(1.5, 0.51, 'TM ≥ 0.50  same fold', fontsize=10, color='#555', va='bottom')
        ax.text(.5, -.17, note, transform=ax.transAxes, ha='center', fontsize=11,
                color='#666', style='italic')
        ax.tick_params(labelsize=11)
    fig.suptitle('Structural resemblance to public TREM2 binders was never screened',
                 fontsize=20, weight='bold', y=1.02)
    fig.text(.5, .965, 'the gate tested sequence identity against the known-binder corpus; '
             'its TM-score criterion applied only to target and control chains',
             ha='center', fontsize=12.5, color='#444')
    plt.tight_layout()
    plt.savefig('fig_novelty_audit_hack.png', dpi=150, facecolor='white'); plt.close()
    print('fig_novelty_audit_hack.png')


# ------------------------------------------------------------------ figure 2
def epitope_map():
    ah = json.load(open('anth_epitopes.json'))
    mu = json.load(open('muni_epitopes.json'))
    hk = json.load(open('epitopes_hack.json'))
    sets = {'Hackathon': hk, 'Muni': mu, 'Anthropic': ah}
    n = {k: len(v) for k, v in sets.items()}
    freq = {k: collections.Counter(r for v in s.values() for r in v) for k, s in sets.items()}

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(15.5, 8.6),
                                  gridspec_kw={'height_ratios': [2.2, 1]})
    xs = sorted(set(itertools.chain.from_iterable(freq[k] for k in freq)))
    xs = [x for x in xs if 15 <= x <= 110]
    w = 0.28
    for i, k in enumerate(['Hackathon', 'Muni', 'Anthropic']):
        ax.bar([x + (i - 1) * w for x in xs], [100 * freq[k].get(x, 0) / n[k] for x in xs],
               width=w, color=COL[k], label=f'{k}  (n={n[k]})')
    ax.set_xlabel('TREM2 residue (156-construct numbering)', fontsize=13)
    ax.set_ylabel('% of designs contacting', fontsize=13)
    ax.set_xlim(15, 110); ax.legend(fontsize=12, frameon=False)
    ax.set_title('All three campaigns converge on the same TREM2 surface',
                 fontsize=16, weight='bold')
    ax.tick_params(labelsize=11)

    # pairwise epitope overlap
    def jac(a, b):
        a, b = set(a), set(b)
        return len(a & b) / max(len(a | b), 1)
    rows = []
    for (ka, sa), (kb, sb) in itertools.combinations(sets.items(), 2):
        v = [jac(x, y) for x in sa.values() for y in sb.values()]
        rows.append((f'{ka} vs {kb}', np.median(v), len(v)))
    for k, s in sets.items():
        v = [jac(x, y) for x, y in itertools.combinations(s.values(), 2)]
        rows.append((f'{k} internal', np.median(v), len(v)))
    rows.sort(key=lambda r: -r[1])
    ax2.barh([r[0] for r in rows], [r[1] for r in rows],
             color=['#888' if 'vs' in r[0] else COL[r[0].split()[0]] for r in rows])
    for i, r in enumerate(rows):
        ax2.text(r[1] + .01, i, f'{r[1]:.2f}', va='center', fontsize=11)
    ax2.set_xlabel('median epitope overlap (Jaccard)', fontsize=13)
    ax2.set_xlim(0, max(r[1] for r in rows) * 1.25)
    ax2.tick_params(labelsize=11)
    plt.tight_layout()
    plt.savefig('fig_epitope_hack.png', dpi=150, facecolor='white'); plt.close()
    print('fig_epitope_hack.png')
    for r in rows:
        print(f'   {r[0]:28s} median Jaccard {r[1]:.3f}  (n={r[2]})')
    top = {k: [r for r, _ in freq[k].most_common(8)] for k in freq}
    print('\n  top-8 contacted residues per collection:')
    for k, v in top.items():
        print(f'   {k:11s} {sorted(v)}')


if __name__ == '__main__':
    novelty_audit()
    epitope_map()

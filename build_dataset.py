#!/usr/bin/env python3
"""
Build every shared artifact the figure scripts consume.

This closes the gap between parse_*.py / pairwise.py and fig_network.py:
  all.fasta, all_binder_pdb/, order_all.json, tm_all.npy, sid_all.npy,
  anth_full.json, archetypes.csv

Run order:  parse_anth.py -> parse_muni.py -> build_dataset.py --stage seqs
            -> (foldseek, see run_pipeline.sh) -> build_dataset.py --stage matrices
"""
import argparse, glob, json, os
import numpy as np, pandas as pd

UP = os.environ.get('UP', '/mnt/user-data/uploads')


# ---------------------------------------------------------------- stage: seqs
def build_fastas():
    m = pd.read_json('muni.json'); a = pd.read_json('anthropic_seqs.json')
    with open('all.fasta', 'w') as f, open('muni.fasta', 'w') as fm, open('anth.fasta', 'w') as fa:
        for _, r in m.iterrows():
            h = f">MUNI|{r['name']}|{r['id']}|{r.get('binding_strength')}"
            f.write(h + "\n" + r.sequence + "\n"); fm.write(h + "\n" + r.sequence + "\n")
        for _, r in a.iterrows():
            h = f">ANTH|{r['full_name']}|{r['workflow']}"
            f.write(h + "\n" + r.sequence + "\n"); fa.write(h + "\n" + r.sequence + "\n")
    print(f'  all.fasta: {len(m) + len(a)} sequences')


def build_pdbs():
    """Binder chains from the Boltz2 co-folds. Chains matched BY SEQUENCE, never by ID:
    the binder is chain B in predicted_boltz2_1to1.cif but chain A in designed.cif."""
    import gemmi
    os.makedirs('all_binder_pdb', exist_ok=True)
    want = {}
    for _, r in pd.read_json('anthropic_seqs.json').iterrows():
        want[r.full_name] = (f'cifs_b2/{r.full_name}/insilico/predicted_boltz2_1to1.cif',
                             r.sequence, r.full_name)
    for _, r in pd.read_json('muni.json').iterrows():
        p = f'muni_files/{r.id}__boltz2_structure_prediction.cif'
        want[r.id] = (p, r.sequence, 'MUNI_' + r.id)
    n = 0
    for key, (path, seq, out) in want.items():
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        st = gemmi.read_structure(path); st.setup_entities()
        chains = {c.name: gemmi.one_letter_code([x.name for x in c]).upper() for c in st[0]}
        hit = [c for c, s in chains.items() if s == seq.upper()]
        if len(hit) != 1:
            raise ValueError(f'{key}: expected exactly one chain matching the design '
                             f'sequence, found {len(hit)}')
        st2 = gemmi.read_structure(path); st2.setup_entities()
        for c in list(st2[0]):
            if c.name != hit[0]:
                st2[0].remove_chain(c.name)
        st2.setup_entities(); st2.write_pdb(f'all_binder_pdb/{out}.pdb'); n += 1
    print(f'  all_binder_pdb/: {n} binder chains')


# ------------------------------------------------------------ stage: matrices
def canonical_order():
    m = pd.read_json('muni.json'); a = pd.read_json('anthropic_seqs.json')
    return ['MUNI_' + x for x in sorted(m.id)] + sorted(a.full_name)


def build_tm(order):
    aln = pd.read_csv('fs_all/aln.tsv', sep='\t', header=None,
                      names=['q', 't', 'tm', 'lddt', 'fident'])
    if len(aln) != len(order) ** 2:
        raise ValueError(f'expected a complete all-vs-all ({len(order)**2} rows), '
                         f'got {len(aln)} — was --exhaustive-search set?')
    R = aln.pivot_table(index='q', columns='t', values='tm', aggfunc='max') \
           .reindex(index=order, columns=order).values
    TM = np.maximum(R, R.T)          # normalise by the shorter chain; see METHODS §3
    np.fill_diagonal(TM, 1.0)
    np.save('tm_all.npy', TM)
    print(f'  tm_all.npy {TM.shape}  symmetric={np.allclose(TM, TM.T)}')
    return TM


def build_sid(order):
    names = json.load(open('names.json')); gid = np.load('gid.npy')
    posM = {n.split('|')[2]: i for i, n in enumerate(names) if n.startswith('MUNI')}
    posA = {n.split('|')[1]: i for i, n in enumerate(names) if n.startswith('ANTH')}
    si = [posM[n[5:]] if n.startswith('MUNI_') else posA[n] for n in order]
    SID = gid[np.ix_(si, si)]
    np.save('sid_all.npy', SID)
    print(f'  sid_all.npy {SID.shape}  symmetric={np.allclose(SID, SID.T)}')
    return SID


def build_anth_full():
    a = pd.read_json('anthropic_seqs.json')
    ds = pd.read_csv(f'{UP}/design_summary.csv', low_memory=False)
    ds = ds[ds.target == 'TREM2']
    cols = ['full_name', 'binder_final', 'kd_nM_final', 'adaptyv_kd_nM', 'twist_kd_nM',
            'vendor_agreement', 'adaptyv_binding', 'twist_binding', 'mouse_binding_final',
            'cyno_binding_final', 'epitope_residues', 'epitope_n_residues',
            'ipsae_min_boltz2', 'sc_dockq_boltz2', 'design_model_status',
            'adaptyv_expression']
    a = a.merge(ds[cols], on='full_name', how='left')
    if a.binder_final.isna().any():
        raise ValueError('binder_final missing for some designs')
    a.to_json('anth_full.json', orient='records')
    print(f'  anth_full.json: {len(a)} rows, {int(a.binder_final.sum())} binders')
    return a


def build_archetypes(order, TM):
    """Fold archetypes: average linkage on 1-TM, cut at 0.40 (TM >= 0.60)."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    a = pd.read_json('anth_full.json').set_index('full_name')
    mu = pd.read_json('muni.json').set_index('id')
    D = 1 - TM.copy(); np.fill_diagonal(D, 0); D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), 'average')
    c = fcluster(Z, t=0.40, criterion='distance')
    df = pd.DataFrame(dict(design=order, cluster=c))
    df['coll'] = np.where(df.design.str.startswith('MUNI_'), 'Muni', 'Anthropic')
    df['workflow'] = df.design.map(lambda x: 'muni' if x.startswith('MUNI_')
                                   else a.loc[x, 'workflow'])
    df['binder'] = df.design.map(lambda x: np.nan if x.startswith('MUNI_')
                                 else bool(a.loc[x, 'binder_final']))
    df['kd'] = df.design.map(lambda x: np.nan if x.startswith('MUNI_')
                             else a.loc[x, 'kd_nM_final'])
    sz = df.cluster.value_counts()
    mapping = {o: f'F{i + 1}' for i, o in enumerate(sz[sz >= 5].index)}
    df['arch'] = df.cluster.map(lambda x: mapping.get(x, 'singletons'))
    df.to_csv('archetypes.csv', index=False)
    print(f'  archetypes.csv: {df.arch.nunique()} archetype labels')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['seqs', 'matrices', 'all'], default='all')
    args = ap.parse_args()
    if args.stage in ('seqs', 'all'):
        print('stage: seqs'); build_fastas(); build_pdbs()
    if args.stage in ('matrices', 'all'):
        print('stage: matrices')
        order = canonical_order()
        json.dump(order, open('order_all.json', 'w'))
        TM = build_tm(order); build_sid(order); build_anth_full()
        build_archetypes(order, TM)

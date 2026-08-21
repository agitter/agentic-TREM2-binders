#!/usr/bin/env python3
"""
Three-collection dataset: Adaptyv x MUNI hackathon (100) + Muni autoresearch (10)
+ Anthropic TREM2 (90) = 200 designs.

Artifacts are suffixed _hack so nothing from the two-collection analysis is touched.

TREM2 construct note: the hackathon used TWO constructs, a 126-mer (92 designs) and
the 156-mer (8 designs).  The 126-mer carries a 2-residue TG extension at the N
terminus, so residue N in 126-numbering == residue N-2 in 156-numbering.  Muni and
Anthropic's long construct use the 156-mer.  Epitope residues are mapped onto
156-numbering here, and every mapping is verified against the residue identity in
the structure rather than trusted from the offset alone.
"""
import argparse, csv, glob, json, os, collections
import numpy as np, pandas as pd, gemmi

csv.field_size_limit(10 ** 9)
UP = os.environ.get('UP', '/mnt/user-data/uploads')
HACK_CSV = os.path.join(UP, 'proteinbase_collection_adaptyv-x-muni-hackathon-ai-agents-vs-humans.csv')
HACK_CIFS = 'hack_cifs'
BINDER_LABELS = {'Strong', 'Medium', 'Weak'}   # == the `binding` boolean, verified 100/100


def read_hack():
    rows = []
    with open(HACK_CSV, newline='', encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            ev = json.loads(r['evaluations'])
            agg = collections.defaultdict(list)
            for e in ev:
                agg[e['metric']].append(e.get('value'))
            first = lambda k: (agg[k][0] if agg.get(k) else None)
            rec = dict(design='HACK_' + r['id'], pb_id=r['id'], name=r['name'],
                       team=r['author'], design_method=r.get('designMethod') or None,
                       sequence=r['sequence'].upper(), length=len(r['sequence']),
                       binding_strength=first('binding_strength'),
                       binding=first('binding'), expressed=first('expressed'),
                       design_class=first('design_class'),
                       classification=first('classification'),
                       kd_M=first('kd'))
            rec['binder'] = rec['binding_strength'] in BINDER_LABELS
            rec['interface'] = next((e['value'] for e in ev
                                     if e['metric'] == 'interface_residues'), None)
            rows.append(rec)
    return pd.DataFrame(rows)


def extract_binder_pdbs(hack):
    """Binder chains for all three collections into one directory (match by sequence)."""
    os.makedirs('all_binder_pdb_hack', exist_ok=True)
    for f in glob.glob('all_binder_pdb/*.pdb'):        # Muni + Anthropic, already verified
        dst = 'all_binder_pdb_hack/' + os.path.basename(f)
        if not os.path.exists(dst):
            open(dst, 'wb').write(open(f, 'rb').read())
    tgt_len = {}
    for _, r in hack.iterrows():
        p = f'{HACK_CIFS}/{r.pb_id}__boltz2_structure_prediction.cif'
        st = gemmi.read_structure(p); st.setup_entities()
        ch = {c.name: gemmi.one_letter_code([x.name for x in c]).upper() for c in st[0]}
        hit = [c for c, s in ch.items() if s == r.sequence]
        if len(hit) != 1:
            raise ValueError(f'{r.pb_id}: {len(hit)} chains match the design sequence')
        other = [c for c in ch if c != hit[0]][0]
        tgt_len[r.design] = len(ch[other])
        st2 = gemmi.read_structure(p); st2.setup_entities()
        for c in list(st2[0]):
            if c.name != hit[0]:
                st2[0].remove_chain(c.name)
        st2.setup_entities()
        st2.write_pdb(f'all_binder_pdb_hack/{r.design}.pdb')
    return tgt_len


def target_seq(pb_id, binder_seq):
    st = gemmi.read_structure(f'{HACK_CIFS}/{pb_id}__boltz2_structure_prediction.cif')
    st.setup_entities()
    ch = {c.name: c for c in st[0]}
    b = [c for c, cc in ch.items()
         if gemmi.one_letter_code([x.name for x in cc]).upper() == binder_seq][0]
    t = [c for c in ch if c != b][0]
    return {r.seqid.num: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
            for r in ch[t]}


def hack_epitopes(hack, tgt_len):
    """Target-side epitope residues mapped onto 156-numbering, identity-verified.

    The 126-mer is TG + TREM2(1..115) + GTKHHHHHH: a 2-residue N-terminal extension
    and a linker+His6 purification tag.  Real TREM2 positions map with a -2 shift;
    contacts to the tag (126-numbering >= 118) are dropped, since they are contacts
    to the construct rather than to the antigen.
    """
    TAG_START_126 = 118
    ref = None
    for _, r in hack.iterrows():
        if tgt_len[r.design] == 156:
            ref = target_seq(r.pb_id, r.sequence); break
    ep, bad, dropped = {}, 0, 0
    for _, r in hack.iterrows():
        if not r.interface:
            continue
        is126 = tgt_len[r.design] == 126
        off = 2 if is126 else 0
        keep = []
        local = target_seq(r.pb_id, r.sequence)
        for x in r.interface:
            if x['chain'] != 'A':
                continue
            n = int(x['residue'])
            if is126 and n >= TAG_START_126:
                dropped += 1
                continue
            m = n - off
            if ref.get(m) and local.get(n) and ref[m] != local[n]:
                bad += 1
            keep.append(m)
        ep[r.design] = sorted(set(keep))
    print(f'  epitope residues remapped to 156-numbering; '
          f'tag contacts dropped: {dropped}; identity mismatches: {bad}')
    return ep


def build_fasta_and_order(hack):
    m = pd.read_json('muni.json'); a = pd.read_json('anthropic_seqs.json')
    order = (sorted('HACK_' + hack.pb_id) + ['MUNI_' + x for x in sorted(m.id)]
             + sorted(a.full_name))
    seqs = {}
    for _, r in hack.iterrows(): seqs[r.design] = r.sequence
    for _, r in m.iterrows(): seqs['MUNI_' + r.id] = r.sequence
    for _, r in a.iterrows(): seqs[r.full_name] = r.sequence
    with open('all_hack.fasta', 'w') as f:
        for n in order:
            f.write(f'>{n}\n{seqs[n]}\n')
    json.dump(order, open('order_hack.json', 'w'))
    json.dump(seqs, open('seqs_hack.json', 'w'))
    return order, seqs


def build_tm(order):
    aln = pd.read_csv('fs_hack/aln.tsv', sep='\t', header=None,
                      names=['q', 't', 'tm', 'lddt', 'fident'])
    if len(aln) != len(order) ** 2:
        raise ValueError(f'expected {len(order)**2} rows, got {len(aln)}')
    R = aln.pivot_table(index='q', columns='t', values='tm', aggfunc='max') \
           .reindex(index=order, columns=order).values
    TM = np.maximum(R, R.T); np.fill_diagonal(TM, 1.0)
    np.save('tm_hack.npy', TM)
    print(f'  tm_hack.npy {TM.shape} symmetric={np.allclose(TM, TM.T)}')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--stage', default='all',
                                                    choices=['seqs', 'matrices', 'all'])
    args = ap.parse_args()
    hack = read_hack()
    print(f'hackathon: {len(hack)} designs, {int(hack.binder.sum())} binders '
          f'(Strong/Medium/Weak), {hack.design_class.nunique()} design classes')
    if args.stage in ('seqs', 'all'):
        tgt_len = extract_binder_pdbs(hack)
        print('  target construct:', dict(collections.Counter(tgt_len.values())))
        hack['target_len'] = hack.design.map(tgt_len)
        json.dump(hack_epitopes(hack, tgt_len), open('epitopes_hack.json', 'w'))
        order, _ = build_fasta_and_order(hack)
        print(f'  all_hack.fasta: {len(order)} sequences; '
              f'{len(glob.glob("all_binder_pdb_hack/*.pdb"))} binder PDBs')
        hack.drop(columns=['interface']).to_csv('hackathon_designs.csv', index=False)
    if args.stage in ('matrices', 'all'):
        build_tm(json.load(open('order_hack.json')))


if __name__ == '__main__':
    main()

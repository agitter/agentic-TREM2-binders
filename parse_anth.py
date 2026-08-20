import re, json, pandas as pd
import os
UP = os.environ.get('UP', '.')
recs=[]; name=None; seq=[]
def flush():
    if name: recs.append((name, ''.join(seq)))
for line in open(os.path.join(UP, 'designs.fasta')):
    line=line.rstrip()
    if line.startswith('>'):
        flush(); name=line[1:]; seq=[]
    else: seq.append(line)
flush()
rows=[]
for hdr,s in recs:
    parts=hdr.split()
    kv=dict(p.split('=',1) for p in parts[1:] if '=' in p)
    if kv.get('target','').upper()!='TREM2': continue
    fid=parts[0]
    m=re.match(r'(.+)_rank(\d+)$', fid)
    stem, rank = m.group(1), int(m.group(2))
    # workflow = model arm + campaign
    wm=re.match(r'(mythos_preview|opus_4_8)_(multi|single)_target_trem2', stem)
    rows.append(dict(full_name=fid, uuid=kv.get('uuid'), stem=stem, rank=rank,
                     model=wm.group(1), campaign=wm.group(2)+'-target',
                     workflow=wm.group(1)+'_'+wm.group(2),
                     sequence=s, length=len(s), hdr_length=int(kv.get('length',0))))
df=pd.DataFrame(rows)
assert (df.length==df.hdr_length).all()
print(df.groupby(['model','campaign']).agg(n=('full_name','size'), len_mean=('length','mean'),
      len_min=('length','min'), len_max=('length','max')))
df.to_json('anthropic_seqs.json',orient='records',indent=1)
print('\ntotal', len(df), 'unique seqs', df.sequence.nunique())

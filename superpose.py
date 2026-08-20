import numpy as np, pandas as pd, re, glob, os, json
def ca_coords(p):
    xs=[]
    for l in open(p):
        if l.startswith('ATOM') and l[12:16].strip()=='CA':
            xs.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return np.array(xs)
CA={os.path.basename(f)[:-4]:ca_coords(f) for f in glob.glob('all_binder_pdb/*.pdb')}
aln=pd.read_csv('fs_sup/sup.tsv',sep='\t',header=None,
    names=['q','t','tm','qs','qe','ts','te','cigar'])
def pairs(qs,ts,cigar):
    qi,ti=qs-1,ts-1; P=[]
    for n,op in re.findall(r'(\d+)([MID])',cigar):
        n=int(n)
        if op=='M':
            P+= [(qi+k,ti+k) for k in range(n)]; qi+=n; ti+=n
        elif op=='I': qi+=n
        else: ti+=n
    return P
def kabsch(P,Q):
    pc,qc=P.mean(0),Q.mean(0)
    H=(P-pc).T@(Q-qc); U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T))
    R=Vt.T@np.diag([1,1,d])@U.T
    return R,qc,pc
def superpose_onto(ref,mob,row):
    P=pairs(int(row.qs),int(row.ts),row.cigar)
    if len(P)<4: return None
    A=CA[ref][[i for i,_ in P]]; B=CA[mob][[j for _,j in P]]
    if len(A)!=len(B): return None
    R,qc,pc=kabsch(B,A)   # rotate mob(B) onto ref(A)
    return (CA[mob]-B.mean(0))@R.T + A.mean(0)
np.save('ca_lens.npy',np.array([len(v) for v in CA.values()]))
json.dump({k:v.tolist() for k,v in CA.items()},open('ca.json','w'))
aln.to_json('aln_sup.json')
print('CA loaded for',len(CA),'structures; lengths',min(len(v) for v in CA.values()),'-',max(len(v) for v in CA.values()))

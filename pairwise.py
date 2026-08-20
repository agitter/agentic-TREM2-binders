import numpy as np, pandas as pd, itertools, json
from Bio import Align
from Bio.Align import substitution_matrices

names=[]; seqs=[]
for line in open('all.fasta'):
    line=line.strip()
    if line.startswith('>'): names.append(line[1:])
    else: seqs.append(line)
n=len(names)

blosum=substitution_matrices.load("BLOSUM62")
# Global (Needleman-Wunsch) aligner, standard EMBOSS-like gap costs
gl=Align.PairwiseAligner(); gl.mode='global'; gl.substitution_matrix=blosum
gl.open_gap_score=-11; gl.extend_gap_score=-1
gl.target_end_gap_score=0; gl.query_end_gap_score=0   # semi-global: free end gaps
# Local (Smith-Waterman)
lo=Align.PairwiseAligner(); lo.mode='local'; lo.substitution_matrix=blosum
lo.open_gap_score=-11; lo.extend_gap_score=-1

def ident(al):
    a,b = al[0], al[1]
    match=sum(1 for x,y in zip(a,b) if x==y and x!='-')
    alnlen=len(a)
    ungapped=sum(1 for x,y in zip(a,b) if x!='-' and y!='-')
    return match, alnlen, ungapped

gid=np.zeros((n,n)); lid=np.zeros((n,n)); lsc=np.zeros((n,n)); lcov=np.zeros((n,n))
for i in range(n):
    gid[i,i]=lid[i,i]=100.0; lsc[i,i]=lo.score(seqs[i],seqs[i]); lcov[i,i]=1.0
for i,j in itertools.combinations(range(n),2):
    a=gl.align(seqs[i],seqs[j])[0]
    m,alen,ung=ident(a)
    # identity normalised by shorter sequence (standard for variable-length binders)
    v=100.0*m/min(len(seqs[i]),len(seqs[j]))
    gid[i,j]=gid[j,i]=v
    b=lo.align(seqs[i],seqs[j])[0]
    m2,alen2,ung2=ident(b)
    lid[i,j]=lid[j,i]=100.0*m2/max(alen2,1)
    lsc[i,j]=lsc[j,i]=b.score
    lcov[i,j]=lcov[j,i]=alen2/min(len(seqs[i]),len(seqs[j]))
np.save('gid.npy',gid); np.save('lid.npy',lid); np.save('lsc.npy',lsc); np.save('lcov.npy',lcov)
json.dump(names,open('names.json','w'))
json.dump(seqs,open('seqs.json','w'))
print('done', gid.shape)

import pandas as pd, json, collections
import os
UP = os.environ.get('UP', '.')
df = pd.read_csv(os.path.join(UP, 'proteinbase_collection_muni-proteina-complex-auto-research.csv'),
                 encoding='utf-8-sig')
rows=[]
for _,r in df.iterrows():
    ev = json.loads(r['evaluations'])
    rec = {'id':r['id'],'name':r['name'],'sequence':r['sequence'],'author':r['author'],
           'length':len(r['sequence'])}
    agg=collections.defaultdict(list)
    for e in ev:
        agg[(e['metric'], e.get('type'))].append(e.get('value'))
    for (m,t),vals in agg.items():
        if m in ('binding_strength','binding','kd','kon','koff','expressed','design_class',
                 'foldstring','classification','novelty','ted_confidence',
                 'boltz2_iptm','boltz2_ipsae','boltz2_plddt','esmfold_plddt',
                 'boltz2_pdockq','boltz2_lis','proteinmpnn_score','molecular_weight',
                 'isoelectric_point','boltz2_complex_iplddt','boltz2_ptm',
                 'shape_complimentarity_boltz2_binder_ss','boltz2_min_ipsae'):
            rec[m+'_all']=vals
            # representative
            if isinstance(vals[0],(int,float)) and not isinstance(vals[0],bool):
                rec[m]=sum(vals)/len(vals); rec[m+'_n']=len(vals)
            else:
                rec[m]=vals[0]; rec[m+'_n']=len(vals)
        if m=='boltz2_structure_prediction':
            rec['struct_url']=vals[0]['url']
        if m=='interface_residues':
            rec['n_interface_res']=len(vals[0])
    rows.append(rec)
m = pd.DataFrame(rows)
m['collection']='Muni'
m['workflow']='muni_proteina'
m.to_json('muni.json',orient='records',indent=1)
print(m[['id','name','length','binding_strength','binding_strength_n','kd','binding','expressed','foldstring','design_class','boltz2_iptm','esmfold_plddt']].to_string())
print()
print(m.binding_strength.value_counts())

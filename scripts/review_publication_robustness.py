from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from run_analysis import ROOT,P,R,design,ols_all,sensitivities,dump

genes=json.loads((R/'RNAseq_candidate_lock.json').read_text())['genes']
discovery=pd.read_csv(R/'exploratory_three_cohort_intersection.tsv',sep='\t').set_index('gene').discovery_beta
rows=[]
for label in ['GSE63514','GSE138080','GSE7803']:
 c_all=pd.read_csv(P/f'{label}_metadata.tsv',sep='\t',index_col=0);c=c_all.loc[c_all.primary_included]
 e=pd.read_csv(P/f'{label}_gene_log_expression.tsv.gz',sep='\t',index_col=0)
 s=pd.read_csv(P/f'{label}_proliferation_scores.tsv',sep='\t',index_col=0)
 rows+=sensitivities(label,c,e,s,genes)
 ranks=e.loc[genes,c.index].rank(axis=1,pct=True)
 rr=ols_all(ranks,design(c,s))
 for row in rr.to_dict('records'):row.update({'cohort':label,'sensitivity':'rank_expression_outcome','n':len(c),'cancer':int(c.group.sum())});rows.append(row)
 if label=='GSE63514':
  for name,mask in [('all_sections_with_dissection_adjustment',c_all.stage.isin(['CIN3','Cancer'])),('CIN2_CIN3_laser_only',c_all.stage.isin(['CIN2','CIN3','Cancer'])&c_all.whole_section.eq(0))]:
   cc=c_all.loc[mask];rr=ols_all(e.loc[genes],design(cc,s,dissection=name.startswith('all_sections')))
   for row in rr.to_dict('records'):row.update({'cohort':label,'sensitivity':name,'n':len(cc),'cancer':int(cc.group.sum())});rows.append(row)
d=pd.DataFrame(rows);d['direction_concordant']=[np.sign(b)==np.sign(discovery[g]) for g,b in zip(d.gene,d.beta)]
d.to_csv(R/'three_gene_publication_robustness_all.tsv',sep='\t',index=False)
summary=[]
for (cohort,gene),g in d.groupby(['cohort','gene']):
 loo=g[g.sensitivity.str.startswith('leave_out_')];other=g[~g.sensitivity.str.startswith('leave_out_')]
 summary.append({'cohort':cohort,'gene':gene,'LOO_models':len(loo),'LOO_direction_concordant_fraction':loo.direction_concordant.mean(),'LOO_beta_min':loo.beta.min(),'LOO_beta_max':loo.beta.max(),'LOO_nominal_p_lt_05_fraction':loo.p.lt(.05).mean(),'other_models':len(other),'other_direction_concordant_fraction':other.direction_concordant.mean(),'other_nominal_p_lt_05_count':int(other.p.lt(.05).sum())})
summary=pd.DataFrame(summary);summary.to_csv(R/'three_gene_publication_robustness_summary.tsv',sep='\t',index=False)
print(summary.to_string(index=False));print('\nNon-LOO checks:');print(d[~d.sensitivity.str.startswith('leave_out_')][['cohort','gene','sensitivity','beta','p','direction_concordant']].to_string(index=False))

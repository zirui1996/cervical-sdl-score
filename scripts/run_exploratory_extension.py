from pathlib import Path
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,P,R,design,ols_all,dump

d=pd.read_csv(R/'discovery_all_genes.tsv',sep='\t');selected=d[d.fdr.lt(.05)&d.beta.abs().ge(.5)].copy();genes=selected.gene.tolist()
rows=[]
for label in ['GSE138080','GSE7803']:
 c=pd.read_csv(P/f'{label}_metadata.tsv',sep='\t',index_col=0);c=c.loc[c.primary_included]
 e=pd.read_csv(P/f'{label}_gene_log_expression.tsv.gz',sep='\t',index_col=0);s=pd.read_csv(P/f'{label}_proliferation_scores.tsv',sep='\t',index_col=0)
 measured=[g for g in genes if g in e.index and e.loc[g,c.index].notna().all()]
 stats=ols_all(e.loc[measured],design(c,s)).set_index('gene')
 for _,candidate in selected.iterrows():
  row={'cohort':label,'gene':candidate.gene,'discovery_beta':candidate.beta,'discovery_fdr':candidate.fdr,'n':len(c)}
  if candidate.gene in stats.index:row.update(stats.loc[candidate.gene].to_dict());row['status']='tested'
  else:row.update({'p':1.,'beta':np.nan,'status':'not_measurable_complete_case'})
  row['direction_concordant']=np.sign(row['beta'])==np.sign(candidate.beta);rows.append(row)
v=pd.DataFrame(rows);v['extension_fdr']=v.groupby('cohort').p.transform(lambda x:multipletests(x,method='fdr_bh')[1]);v['extension_supported']=v.extension_fdr.lt(.05)&v.direction_concordant
v.to_csv(R/'exploratory_all_qualified_external.tsv',sep='\t',index=False)
supported=v.groupby('gene').extension_supported.sum();both=supported[supported.eq(2)].index.tolist()
ranked=[]
for gene in both:
 ss=v[v.gene.eq(gene)];ranked.append({'gene':gene,'worst_external_fdr':ss.extension_fdr.max(),'discovery_p':float(selected.set_index('gene').loc[gene,'p']),'discovery_beta':float(selected.set_index('gene').loc[gene,'beta'])})
if ranked:
 ranked=pd.DataFrame(ranked).sort_values(['worst_external_fdr','discovery_p','gene']);ranked.to_csv(R/'exploratory_three_cohort_intersection.tsv',sep='\t',index=False)
else:ranked=pd.DataFrame(columns=['gene','worst_external_fdr','discovery_p','discovery_beta']);ranked.to_csv(R/'exploratory_three_cohort_intersection.tsv',sep='\t',index=False)
summary={'analysis_status':'post hoc exploratory; original top-five screen remains failed','qualified_discovery_pool':len(genes),'supported_primary_external':int(v.loc[v.cohort.eq('GSE138080'),'extension_supported'].sum()),'supported_secondary_external':int(v.loc[v.cohort.eq('GSE7803'),'extension_supported'].sum()),'intersection_count':len(both),'intersection_genes':ranked.gene.tolist()}
dump('exploratory_extension_summary.json',summary);print(json.dumps(summary,indent=2));print(ranked.head(10).to_string(index=False))

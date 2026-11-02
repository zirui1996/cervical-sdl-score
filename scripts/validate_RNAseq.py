from pathlib import Path
import csv,gzip,json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,RAW,P,R,F,SEED,design,ols_all,scores,gene_sets,diagnostics,sensitivities,auc_bootstrap,dump
import matplotlib.pyplot as plt

lock=json.loads((R/'RNAseq_candidate_lock.json').read_text());genes=lock['genes']
assert genes==json.loads((R/'exploratory_extension_summary.json').read_text())['intersection_genes']
m={}
with gzip.open(RAW/'GSE223804_series_matrix.txt.gz','rt') as h:
 for line in h:
  if line.startswith('!Sample_'):
   v=next(csv.reader([line],delimiter='\t'));m.setdefault(v[0],v[1:])
  if line.startswith('!series_matrix_table_begin'):break
c=pd.DataFrame({'sample':m['!Sample_geo_accession'],'title':m['!Sample_title']})
c['matrix_column']=c.title.str.extract(r'\[([^\]]+)\]',expand=False).str.replace('_','',regex=False)
c['group']=c.matrix_column.str.startswith('CR').astype(int)
c['stage']=np.where(c.group.eq(1),'Cancer',np.where(c.matrix_column.eq('CIN5'),'CIN1',np.where(c.matrix_column.str.startswith('CIN'),'HSIL','Normal')))
c['primary_included']=c.stage.isin(['Cancer','HSIL']);c=c.set_index('matrix_column');c.to_csv(P/'GSE223804_metadata.tsv',sep='\t')
raw=pd.read_csv(RAW/'GSE223804_cpm_CIN_CR_Norm.txt.gz',sep='\t',index_col=0)
assert set(c.index)==set(raw.columns) and raw.values.min()>=0
e=np.log2(raw+1);e=e.groupby(level=0).median()
cache=P/'GSE223804_log2_CPMplus1.tsv.gz';tmp=P/'GSE223804_log2_CPMplus1.atomic.tsv.gz'
e.to_csv(tmp,sep='\t',compression='gzip');reopened=pd.read_csv(tmp,sep='\t',index_col=0)
assert reopened.shape==e.shape and np.allclose(reopened,e)
tmp.replace(cache)
common=pd.read_csv(R/'platform_common_gene_universe.tsv',sep='\t').gene.tolist();s,cov=scores(e,common,gene_sets());s.to_csv(P/'GSE223804_proliferation_scores.tsv',sep='\t');cov.to_csv(R/'GSE223804_score_coverage.tsv',sep='\t',index=False)
cc=c.loc[c.primary_included];assert cc.group.value_counts().to_dict()=={1:9,0:4}
direction=pd.read_csv(R/'exploratory_three_cohort_intersection.tsv',sep='\t').set_index('gene').discovery_beta
res=ols_all(e.loc[genes],design(cc,s));res['holm_p']=multipletests(res.p,method='holm')[1]
res['direction_concordant']=[np.sign(b)==np.sign(direction[g]) for g,b in zip(res.gene,res.beta)];res['replicated']=res.holm_p.lt(.05)&res.direction_concordant
ua=ols_all(e.loc[genes],pd.DataFrame({'const':1.,'group':cc.group},index=cc.index));res=res.merge(ua[['gene','beta','p']],on='gene',suffixes=('','_unadjusted'))
for gene in genes:
 metrics=auc_bootstrap(e.loc[gene,cc.index].values,cc.group.values,np.sign(direction[gene]))
 for k,v in metrics.items():res.loc[res.gene.eq(gene),k]=v
res.to_csv(R/'new_RNAseq_validation.tsv',sep='\t',index=False)
sen=sensitivities('GSE223804',cc,e,s,genes)
inc=c.loc[c.stage.isin(['HSIL','CIN1','Cancer'])]
extra=ols_all(e.loc[genes],design(inc,s))
for row in extra.to_dict('records'):row.update({'cohort':'GSE223804','sensitivity':'includes_low_grade_CIN1','n':len(inc),'cancer':int(inc.group.sum())});sen.append(row)
pd.DataFrame(sen).to_csv(R/'new_RNAseq_sensitivities.tsv',sep='\t',index=False);dump('new_RNAseq_design_diagnostics.json',diagnostics('GSE223804',cc,s))
dump('new_RNAseq_summary.json',{'n':len(cc),'high_grade_precancer':4,'invasive_SCC':9,'excluded_normal':2,'excluded_CIN1':1,'genes':genes,'replicated':res.loc[res.replicated,'gene'].tolist(),'gene_expression_scale':'log2(author provided CPM + 1)','source_file_column_sums':raw.sum().to_dict(),'status':'new holdout validation after explicitly post hoc three-microarray extension'})
fig,axs=plt.subplots(1,len(genes),figsize=(13,4))
for ax,gene in zip(axs,genes):
 for k,col in [(0,'#457b9d'),(1,'#d56c4a')]:
  ix=cc.index[cc.group.eq(k)];ax.scatter(s.loc[ix,'Composite'],e.loc[gene,ix],c=col,label='SCC' if k else 'High-grade CIN',s=35)
 ax.set(title=gene,xlabel='Proliferation rank score',ylabel='log2(author CPM + 1)');ax.legend(fontsize=8)
fig.suptitle('New RNA-seq holdout: 4 high-grade precancers versus 9 invasive SCC');fig.tight_layout();fig.savefig(F/'05_new_RNAseq_holdout.png',dpi=150);fig.savefig(F/'05_new_RNAseq_holdout.pdf');plt.close(fig)
print(res.to_string(index=False));print('PASSED',res.loc[res.replicated,'gene'].tolist())

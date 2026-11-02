"""Post hoc revision checks (2026-10-30); the frozen ten-gene score is not changed.

Uses only saved processed matrices and saved proliferation scores, because the
raw GMT/series files are not part of the manuscript package. The fixed score is
recomputed from the lock and must reproduce the published adjusted estimates
before any new output is written.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data/processed';R=ROOT/'results';T=ROOT/'manuscript/tables'
LOCK=json.loads((R/'joint_program_lock.json').read_text())
UNIVERSE=pd.read_csv(R/'platform_common_gene_universe.tsv',sep='\t').gene.tolist()
DIFF=['KRT4','KRT13','IVL','SPRR1A','SPRR3','KRT1','DSG1']
BASAL=['KRT5','KRT14','TP63','DSG3','DSC3','SFN']
GLAND=['PAX8','TFF3','AGR2','MUC5B','MUC5AC','KRT7']

def ols_hc3(y,X):
 A=X.values.astype(float);Y=y.values.astype(float);inv=np.linalg.inv(A.T@A);b=inv@A.T@Y;res=Y-A@b
 h=np.einsum('ij,jk,ik->i',A,inv,A);w=(A@inv)[:,1];se=np.sqrt(np.sum((res/(1-h))**2*w**2));df=len(A)-A.shape[1]
 t=b[1]/se;crit=stats.t.ppf(.975,df);return {'beta':b[1],'se_HC3':se,'p':2*stats.t.sf(abs(t),df),'ci_low':b[1]-crit*se,'ci_high':b[1]+crit*se,'residual_df':df}

def z(v):return (v-v.mean())/v.std(ddof=1)

def ranks(e):return e.loc[e.index.intersection(UNIVERSE)].rank(axis=0,pct=True)

def signed_score(rank,genes,directions,minimum):
 have=[g for g in genes if g in rank.index]
 if len(have)<minimum:return None,len(have)
 return sum(directions[g]*rank.loc[g] for g in have)/len(have),len(have)

def load():
 out={}
 for label in ['GSE63514','GSE138080','GSE7803']:
  c=pd.read_csv(P/f'{label}_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/f'{label}_gene_log_expression.tsv.gz',sep='\t',index_col=0)
  s=pd.read_csv(P/f'{label}_proliferation_scores.tsv',sep='\t',index_col=0);out[label]=(c,e,s,'Development')
 for label,role in [('GSE166466','Additional microarray'),('GSE149763','Additional RNA-seq'),('GSE147009','Additional RNA-seq')]:
  c=pd.read_csv(P/f'{label}_extension_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/f'{label}_extension_gene_log_expression.tsv.gz',sep='\t',index_col=0)
  s=pd.read_csv(P/f'{label}_extension_scores.tsv',sep='\t',index_col=0);out[label]=(c,e,s,role)
 # The packaged GSE14404 gene-intensity gzip is truncated; only its saved
 # frozen/proliferation scores are usable, so top-k and marker checks skip it.
 c=pd.read_csv(P/'GSE14404_independent_metadata.tsv',sep='\t',index_col=0)
 s=pd.read_csv(P/'GSE14404_raw_Cy5_scores.tsv',sep='\t',index_col=0);out['GSE14404']=(c,None,s,'Additional microarray')
 c=pd.read_csv(P/'GSE223804_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/'GSE223804_log2_CPMplus1.tsv.gz',sep='\t',index_col=0)
 s=pd.read_csv(P/'GSE223804_proliferation_scores.tsv',sep='\t',index_col=0);out['GSE223804']=(c,e,s,'Prior RNA context')
 return out

def main():
 cohorts=load();lockdir={x['gene']:x['direction'] for x in LOCK['candidates']};lockgenes=list(lockdir)
 published=pd.read_csv(T/'S03_fixed_score_tissue_results.tsv',sep='\t').set_index('cohort')
 disc=pd.read_csv(R/'meta_discovery_all_genes.tsv',sep='\t')
 q=disc[disc.same_direction_all_three&disc.combined_fdr.lt(.05)&disc.mean_absolute_standardized_beta.ge(.5)].sort_values(['combined_p','gene'])
 assert len(q)==812 and q.gene.head(10).tolist()==lockgenes,'Discovery ordering does not reproduce the lock'
 adj=[];topk=[];corr=[];stage=[]
 for label,(c,e,s,role) in cohorts.items():
  cc=c.loc[c.primary_included.astype(bool)].copy()
  if e is None:rank=None;fixed=s.loc[cc.index,'Joint_program'];m=9
  else:rank=ranks(e)[cc.index];fixed,m=signed_score(rank,lockgenes,lockdir,LOCK['minimum_measured'])
  comp=z(s.loc[cc.index,'Composite']);y=z(fixed);g=cc.group.astype(float)
  Xa=pd.DataFrame({'const':1.,'group':g,'proliferation_z':comp});Xu=Xa[['const','group']]
  a=ols_hc3(y,Xa);u=ols_hc3(y,Xu)
  assert abs(a['beta']-published.loc[label,'beta'])<1e-4,f'{label} fixed score does not reproduce S03'
  if e is None:
   adj.append({'cohort':label,'role':role,'n':len(cc),'measured_members':m,**{'unadjusted_'+k:u[k] for k in ['beta','ci_low','ci_high','p']},**{'adjusted_'+k:a[k] for k in ['beta','ci_low','ci_high','p']},'adjusted_minus_unadjusted':a['beta']-u['beta'],
    **dict(zip(['composite_SCC_minus_HSIL_sd','composite_diff_ci_low','composite_diff_ci_high','composite_diff_p'],[ols_hc3(comp,Xu)[k] for k in ['beta','ci_low','ci_high','p']])),'group_composite_r':stats.pointbiserialr(g,s.loc[cc.index,'Composite']).statistic})
   adj[-1]['group_composite_VIF']=1/(1-adj[-1]['group_composite_r']**2);continue
  r=stats.pointbiserialr(g,s.loc[cc.index,'Composite']).statistic
  pc=ols_hc3(comp,Xu)  # proliferation difference between groups, in Composite SD
  adj.append({'cohort':label,'role':role,'n':len(cc),'measured_members':m,'unadjusted_beta':u['beta'],'unadjusted_ci_low':u['ci_low'],'unadjusted_ci_high':u['ci_high'],'unadjusted_p':u['p'],
   'adjusted_beta':a['beta'],'adjusted_ci_low':a['ci_low'],'adjusted_ci_high':a['ci_high'],'adjusted_p':a['p'],'adjusted_minus_unadjusted':a['beta']-u['beta'],
   'composite_SCC_minus_HSIL_sd':pc['beta'],'composite_diff_ci_low':pc['ci_low'],'composite_diff_ci_high':pc['ci_high'],'composite_diff_p':pc['p'],'group_composite_r':r,'group_composite_VIF':1/(1-r*r)})
  for k in [5,10,15,20]:
   genes=q.gene.head(k).tolist();dirs=dict(zip(q.gene,q.direction.astype(int)))
   v,mm=signed_score(rank,genes,dirs,int(np.ceil(.8*k)))
   rec={'cohort':label,'role':role,'k':k,'members_defined':k,'members_measured':mm,'positive_members':int(sum(dirs[x]>0 for x in genes))}
   if v is not None:rec.update(ols_hc3(z(v),Xa))
   topk.append(rec)
  for name,genes in [('Squamous_differentiated',DIFF),('Squamous_basal',BASAL),('Glandular_secretory',GLAND)]:
   have=rank.index.intersection(genes)
   if len(have)<3:continue
   pv=rank.loc[have].mean(axis=0);rho,p=stats.spearmanr(fixed,pv)
   # Pooled rho can be carried by the HSIL/SCC difference itself, so also
   # report a group-adjusted partial Spearman and within-group rho (n>=6).
   rf=fixed.rank();rp=pv.rank();G=np.column_stack([np.ones(len(g)),g])
   res=lambda v:v-G@np.linalg.lstsq(G,v,rcond=None)[0]
   pr=np.corrcoef(res(rf.values),res(rp.values))[0,1];dfp=len(g)-3;tp=pr*np.sqrt(dfp/(1-pr*pr))
   rec={'cohort':label,'role':role,'proxy':name,'n':len(cc),'measured_markers':len(have),'spearman_rho':rho,'p':p,'group_adjusted_partial_rho':pr,'partial_p':2*stats.t.sf(abs(tp),dfp)}
   for grp,k in [('HSIL',0),('SCC',1)]:
    idx=g.index[g.eq(k)]
    if len(idx)>=6:r2,p2=stats.spearmanr(fixed[idx],pv[idx]);rec.update({f'within_{grp}_n':len(idx),f'within_{grp}_rho':r2,f'within_{grp}_p':p2})
   corr.append(rec)
  # Stage-group marker profile where a normal group exists (descriptive only).
  if 'stage' in c and c.stage.isin(['Normal']).any() and label in ['GSE166466','GSE63514','GSE138080','GSE7803']:
   keep=c.index[c.stage.isin(['Normal','HSIL','CIN3','Cancer'])&c.index.isin(e.columns)]
   if label=='GSE63514':keep=c.index[c.stage.isin(['Normal','CIN3','Cancer'])&c.whole_section.eq(0)]
   if label=='GSE7803':keep=keep[~c.loc[keep,'excluded_test'].astype(bool)]
   rk=ranks(e[keep]);fs,_=signed_score(rk,lockgenes,lockdir,8)
   for st,idx in c.loc[keep].groupby('stage').groups.items():
    row={'cohort':label,'stage':'HSIL' if st=='CIN3' else st,'n':len(idx),'fixed_score_mean_raw':float(fs[idx].mean())}
    for name,genes in [('Squamous_differentiated',DIFF),('Squamous_basal',BASAL),('Glandular_secretory',GLAND)]:
     row[name+'_rank_mean']=float(rk.loc[rk.index.intersection(genes),idx].mean(axis=0).mean())
    stage.append(row)
 adj=pd.DataFrame(adj);topk=pd.DataFrame(topk);corr=pd.DataFrame(corr);stage=pd.DataFrame(stage)
 corr['bh_fdr']=multipletests(corr.p,method='fdr_bh')[1];corr['partial_bh_fdr']=multipletests(corr.partial_p,method='fdr_bh')[1]
 independent_stage_checks()
 for d,name in [(adj,'revision_adjusted_vs_unadjusted.tsv'),(topk,'revision_topk_sensitivity.tsv'),(corr,'revision_score_epithelial_marker_correlation.tsv'),(stage,'revision_stage_marker_profiles.tsv')]:
  d.to_csv(R/name,sep='\t',index=False)
 import scipy,statsmodels
 (R/'revision_posthoc_record.json').write_text(json.dumps({'date':'2026-10-30','timing':'post hoc, after all primary results were known; fixed score unchanged','reproduced_S03_adjusted_betas':True,
  'discovery_ordering_reproduces_lock':True,'environment':{'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},
  'note':'Top-k scores reuse the frozen discovery ordering; k=10 is the frozen score. Correlations are descriptive within HSIL+SCC analysis sets.'},indent=2))
 pd.set_option('display.width',250)
 print(adj.round(3).to_string(index=False));print(topk.round(3).to_string(index=False));print(corr.round(3).to_string(index=False));print(stage.round(3).to_string(index=False))

def independent_stage_checks():
 """Grade steps in sources never used for feature selection (saved frozen scores)."""
 import statsmodels.api as sm
 rows=[]
 specs=[('GSE51993',P/'GSE51993_extension_metadata.tsv',P/'GSE51993_extension_scores.tsv',['Normal','CIN1','HSIL']),
        ('GSE14404',P/'GSE14404_independent_metadata.tsv',P/'GSE14404_raw_Cy5_scores.tsv',['Normal','CIN1','HSIL','Cancer'])]
 for label,mp,sp,levels in specs:
  c=pd.read_csv(mp,sep='\t',index_col=0);s=pd.read_csv(sp,sep='\t',index_col=0)
  cc=c.loc[c.stage.isin(levels)].copy()
  if label=='GSE14404':cc=cc.loc[cc.stage.ne('Cancer')|cc.primary_included.astype(bool)]  # explicitly annotated SCC only
  y=z(s.loc[cc.index,'Joint_program']);X=pd.DataFrame({'const':1.},index=cc.index)
  for lv in levels[1:]:X[lv]=cc.stage.eq(lv).astype(float)
  X['proliferation_z']=z(s.loc[cc.index,'Composite'])
  fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);names=list(X.columns)
  def vec(pos,neg=None):
   v=np.zeros(len(names))
   if pos in names:v[names.index(pos)]+=1
   if neg in names:v[names.index(neg)]-=1
   return v
  steps=[(levels[i],levels[i-1]) for i in range(1,len(levels))]
  contrasts=[(f'{a}_vs_{b}',vec(a,b if b!=levels[0] else None)) for a,b in steps]
  if len(levels)==4:contrasts.append(('Cancer_vs_HSIL_minus_HSIL_vs_Normal',vec('Cancer','HSIL')-vec('HSIL')))
  for name,v in contrasts:
   t=fit.t_test(v);ci=np.asarray(t.conf_int()).ravel()
   rows.append({'cohort':label,'role':'Not used for feature selection','contrast':name,'n':len(cc),**{f'n_{lv}':int(cc.stage.eq(lv).sum()) for lv in levels},
    'beta':float(np.ravel(t.effect)[0]),'ci_low':ci[0],'ci_high':ci[1],'p':float(np.ravel(t.pvalue)[0]),'residual_df':int(fit.df_resid)})
 d=pd.DataFrame(rows);d.to_csv(R/'revision_independent_stage_checks.tsv',sep='\t',index=False)
 print(d.round(3).to_string(index=False))

if __name__=='__main__':main()

"""Dated exploratory comparison of adjustment and cross-study agreement.

No candidate selection or validation rules are changed. Participant bootstraps
are stratified by pathology, and interval estimates are pointwise/descriptive.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import ROOT,P,R,F,design,ols_all,dump

LABELS=['GSE63514','GSE138080','GSE7803']; PAIRS=[(0,1),(0,2),(1,2)]
SEED=20261030; BOOTSTRAPS=1000

def unadjusted(c):
 return pd.DataFrame({'const':1.,'group':c.group.astype(float)},index=c.index)

def standardized_coefficients(Y,group,score):
 # Y is samples by genes. Recompute both SDs in each bootstrap.
 z=(score-score.mean())/score.std(ddof=1)
 xx=np.column_stack([np.ones(len(group)),group,z])
 inv=np.linalg.inv(xx.T@xx)
 adj=(inv@xx.T@Y)[1]
 raw=Y[group==1].mean(axis=0)-Y[group==0].mean(axis=0)
 sd=Y.std(axis=0,ddof=1)
 return raw/sd,adj/sd

def rho(x,y):
 good=np.isfinite(x)&np.isfinite(y)
 return float(spearmanr(x[good],y[good]).statistic)

def main():
 frozen=pd.read_csv(R/'primary_tested_gene_universe.tsv',sep='\t').gene.tolist()
 cohorts={};statsrows=[];summary=[];available=[]
 for label in LABELS:
  allc=pd.read_csv(P/f'{label}_metadata.tsv',sep='\t',index_col=0);c=allc.loc[allc.primary_included].copy()
  e=pd.read_csv(P/f'{label}_gene_log_expression.tsv.gz',sep='\t',index_col=0)
  s=pd.read_csv(P/f'{label}_proliferation_scores.tsv',sep='\t',index_col=0)
  candidates=e.index.intersection(frozen);v=e.loc[candidates,c.index]
  valid=v.index[v.notna().all(axis=1)&v.std(axis=1).gt(0)]
  available.append(set(valid));cohorts[label]=(allc,c,e,s)
  aa=ols_all(e.loc[valid],design(c,s));uu=ols_all(e.loc[valid],unadjusted(c))
  joined=aa.merge(uu,on='gene',suffixes=('_adjusted','_unadjusted'))
  sd=e.loc[valid,c.index].std(axis=1)
  for kind in ['adjusted','unadjusted']:joined['standardized_beta_'+kind]=joined['beta_'+kind]/joined.gene.map(sd)
  joined.insert(0,'cohort',label);statsrows.append(joined)
  A=joined.fdr_adjusted.lt(.05);U=joined.fdr_unadjusted.lt(.05)
  summary.append({'cohort':label,'genes_tested':len(valid),'unadjusted_fdr_lt_05':int(U.sum()),'adjusted_fdr_lt_05':int(A.sum()),'both':int((A&U).sum()),'unadjusted_only':int((U&~A).sum()),'adjusted_only':int((A&~U).sum()),'coefficient_sign_changes':int((np.sign(joined.beta_adjusted)!=np.sign(joined.beta_unadjusted)).sum()),'median_absolute_standardized_coefficient_change':float((joined.standardized_beta_adjusted-joined.standardized_beta_unadjusted).abs().median())})
 pd.concat(statsrows).to_csv(R/'genomewide_unadjusted_adjusted_coefficients.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(summary).to_csv(R/'genomewide_adjustment_summary.tsv',sep='\t',index=False)
 common=sorted(set.intersection(*available));pd.Series(common,name='gene').to_csv(R/'genomewide_common_measured_universe.tsv',sep='\t',index=False)
 arrays=[];point=[]
 for label in LABELS:
  _,c,e,s=cohorts[label];Y=e.loc[common,c.index].T.to_numpy();g=c.group.to_numpy();score=s.loc[c.index,'Composite'].to_numpy()
  arrays.append((Y,g,score));point.append(standardized_coefficients(Y,g,score))
 rows=[]
 for i,j in PAIRS:
  u=rho(point[i][0],point[j][0]);a=rho(point[i][1],point[j][1])
  rows.append({'study_1':LABELS[i],'study_2':LABELS[j],'genes':len(common),'rho_unadjusted':u,'rho_adjusted':a,'rho_difference_adjusted_minus_unadjusted':a-u,'sign_concordance_unadjusted':float((np.sign(point[i][0])==np.sign(point[j][0])).mean()),'sign_concordance_adjusted':float((np.sign(point[i][1])==np.sign(point[j][1])).mean())})
 rng=np.random.default_rng(SEED);boots=[]
 for b in range(BOOTSTRAPS):
  fitted=[]
  for Y,g,score in arrays:
   ix=np.r_[rng.choice(np.flatnonzero(g==0),int((g==0).sum()),replace=True),rng.choice(np.flatnonzero(g==1),int((g==1).sum()),replace=True)]
   fitted.append(standardized_coefficients(Y[ix],g[ix],score[ix]))
  for i,j in PAIRS:
   u=rho(fitted[i][0],fitted[j][0]);a=rho(fitted[i][1],fitted[j][1])
   boots.append({'bootstrap':b,'study_1':LABELS[i],'study_2':LABELS[j],'rho_unadjusted':u,'rho_adjusted':a,'difference':a-u})
  if (b+1)%200==0:print('Participant bootstrap',b+1,'/',BOOTSTRAPS,flush=True)
 boot=pd.DataFrame(boots);boot.to_csv(R/'genomewide_participant_bootstrap.tsv.gz',sep='\t',index=False,compression='gzip')
 for row in rows:
  bb=boot[(boot.study_1==row['study_1'])&(boot.study_2==row['study_2'])]
  row.update({'difference_ci_low':bb.difference.quantile(.025),'difference_ci_high':bb.difference.quantile(.975),'bootstrap_resamples':len(bb)})
 agree=pd.DataFrame(rows);agree.to_csv(R/'genomewide_cross_study_agreement.tsv',sep='\t',index=False)
 # Discovery sampling analyses retain the entire fixed 6685-gene universe.
 allc,c,e,s=cohorts['GSE63514'];allst=allc.loc[allc.stage.isin(['CIN3','Cancer'])]
 sampling=[]
 for name,cc,XX in [('laser_primary_adjusted',c,design(c,s)),('all_sections_no_dissection_covariate',allst,design(allst,s)),('all_sections_dissection_covariate',allst,design(allst,s,dissection=True))]:
  dd=ols_all(e.loc[frozen],XX);dd.insert(0,'model',name);dd['n']=len(cc);sampling.append(dd)
 cancers=allc.loc[allc.stage.eq('Cancer')].copy();cancers['group']=cancers.whole_section
 for name,XX in [('within_cancer_whole_vs_laser_unadjusted',unadjusted(cancers)),('within_cancer_whole_vs_laser_proliferation_adjusted',design(cancers,s))]:
  dd=ols_all(e.loc[frozen],XX);dd.insert(0,'model',name);dd['n']=len(cancers);sampling.append(dd)
 samp=pd.concat(sampling);samp.to_csv(R/'genomewide_sampling_sensitivity.tsv.gz',sep='\t',index=False,compression='gzip')
 ss=samp.groupby('model').agg(genes=('gene','size'),n=('n','first'),fdr_lt_05=('fdr',lambda x:int((x<.05).sum()))).reset_index();ss.to_csv(R/'genomewide_sampling_summary.tsv',sep='\t',index=False)
 # Explicit numerical verification on representatives of every new model.
 import statsmodels.api as sm
 verify=[]
 for label in LABELS:
  _,cc,ee,sc=cohorts[label]
  for kind,XX in [('unadjusted',unadjusted(cc)),('adjusted',design(cc,sc))]:
   for gene in [common[0],common[len(common)//2],common[-1]]:
    rr=ols_all(ee.loc[[gene]],XX).iloc[0];fit=sm.OLS(ee.loc[gene,cc.index],XX).fit(cov_type='HC3',use_t=True)
    assert np.isclose(rr.beta,fit.params['group']) and np.isclose(rr.p,fit.pvalues['group'],rtol=1e-7,atol=1e-12)
    verify.append({'cohort':label,'model':kind,'gene':gene})
 assert len(boot)==BOOTSTRAPS*3 and boot.select_dtypes('number').notna().all().all()
 record={'extension_plan_sha256':hashlib.sha256((ROOT/'GENOME_WIDE_EXTENSION_v1.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seed':SEED,'bootstraps_per_pair':BOOTSTRAPS,'common_measured_genes':len(common),'cohort_summary':summary,'cross_study_agreement':rows,'numerical_checks':verify,'inference':'Post hoc exploratory; pointwise bootstrap intervals conditional on fixed gene universe; no new marker selection.'}
 dump('genomewide_extension_summary.json',record)
 # Publication figures; no significance stars.
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(1,2,figsize=(11.5,4.1))
 si=pd.DataFrame(summary);x=np.arange(3)
 axs[0].bar(x-.18,si.unadjusted_fdr_lt_05,.36,label='Unadjusted',color='#567c99');axs[0].bar(x+.18,si.adjusted_fdr_lt_05,.36,label='Adjusted',color='#b46f4c');axs[0].set(xticks=x,xticklabels=LABELS,ylabel='Genes with within-model BH FDR < 0.05',title='A  Gene-wise associations');axs[0].legend(frameon=False)
 for k,row in agree.iterrows():
  axs[1].plot([0,1],[row.rho_unadjusted,row.rho_adjusted],'-o',label=row.study_1+' / '+row.study_2)
 axs[1].set(xticks=[0,1],xticklabels=['Unadjusted','Adjusted'],ylabel='Spearman correlation of standardized coefficients',title='B  Cross-study agreement');axs[1].legend(frameon=False,fontsize=7)
 fig.tight_layout();fig.savefig(F/'08_genomewide_adjustment.png',dpi=300);fig.savefig(F/'08_genomewide_adjustment.pdf');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7.3,3.7))
 for k,row in agree.iterrows():ax.errorbar(row.rho_difference_adjusted_minus_unadjusted,k,xerr=[[row.rho_difference_adjusted_minus_unadjusted-row.difference_ci_low],[row.difference_ci_high-row.rho_difference_adjusted_minus_unadjusted]],fmt='o',capsize=3,color='#365a72')
 ax.axvline(0,color='grey',ls='--');ax.set(yticks=range(3),yticklabels=[r.study_1+' / '+r.study_2 for _,r in agree.iterrows()],xlabel='Change in Spearman correlation after adjustment\nPointwise 95% participant-bootstrap interval',title='Adjustment did not have a prespecified beneficial direction');fig.tight_layout();fig.savefig(F/'09_agreement_bootstrap.png',dpi=300);fig.savefig(F/'09_agreement_bootstrap.pdf');plt.close(fig)
 print(si.to_string(index=False));print(agree.to_string(index=False));print(ss.to_string(index=False))

if __name__=='__main__':main()

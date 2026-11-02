import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from run_analysis import R,P,F,design,ols_all,dump
from analyze_stage_and_composition import load_all
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def signed_z(Y,X):
 A=np.asarray(X,float);assert np.linalg.matrix_rank(A)==A.shape[1]
 inv=np.linalg.inv(A.T@A);coef=inv@A.T@Y;res=Y-A@coef;lev=np.einsum('ij,jk,ik->i',A,inv,A);w=(A@inv)[:,1]
 if np.max(lev)>=1-1e-10:raise ValueError('unit_leverage_HC3_not_estimable')
 variance=np.sum((res/(1-lev[:,None]))**2*w[:,None]**2,axis=0);se=np.sqrt(variance)
 with np.errstate(divide='ignore',invalid='ignore'):p=2*stats.t.sf(np.abs(coef[1]/se),len(A)-A.shape[1])
 # A resampled gene may be constant even though it was variable in the full
 # frozen universe. It supplies zero signed evidence, rather than invalidating
 # an entire bootstrap across thousands of other genes.
 constant=np.ptp(Y,axis=0)<=1e-12;p[constant]=1
 return np.sign(coef[1])*stats.norm.isf(np.clip(p,1e-300,1)/2)

def main():
 cohorts=load_all();old=pd.read_csv(R/'genomewide_unadjusted_adjusted_coefficients.tsv.gz',sep='\t');meta=pd.read_csv(R/'meta_discovery_all_genes.tsv',sep='\t').set_index('gene');rng=np.random.default_rng(20261031);summary=[];draws=[];full=[]
 for label in ['GSE166466','GSE149763','GSE147009']:
  c,e,s=cohorts[label];cc=c.loc[c.primary_included];pool=e.loc[e.index.intersection(meta.index),cc.index];genes=pool.index[pool.notna().all(axis=1)&pool.std(axis=1).gt(0)].tolist();Y=e.loc[genes,cc.index].T.to_numpy();X=design(cc,s)
  external=signed_z(Y,X);fit=ols_all(e.loc[genes],X).set_index('gene').loc[genes];reference=np.sign(fit.beta)*stats.norm.isf(np.clip(fit.p,1e-300,1)/2);assert np.allclose(external,reference)
  train={}
  for study in ['GSE63514','GSE138080','GSE7803']:
   dd=old[old.cohort.eq(study)].set_index('gene').loc[genes];train[study]=np.sign(dd.beta_adjusted.to_numpy())*stats.norm.isf(np.clip(dd.p_adjusted.to_numpy(),1e-300,1)/2)
  train['Joint']=meta.loc[genes,'combined_z'].to_numpy();corr={study:stats.spearmanr(v,external).statistic for study,v in train.items()};row={'cohort':label,'n':len(cc),'genes':len(genes),**{'rho_'+k:v for k,v in corr.items()},'delta_joint_minus_single':corr['Joint']-corr['GSE63514']}
  table=pd.DataFrame({'gene':genes,'external_signed_z':external,'external_beta':fit.beta.values,'external_p':fit.p.values,'external_fdr':multipletests(fit.p,method='fdr_bh')[1],'cohort':label});full.append(table)
  group=cc.group.to_numpy();score=s.loc[cc.index,'Composite'].to_numpy();valid=[]
  for b in range(1000):
   ix=np.r_[rng.choice(np.flatnonzero(group==0),int((group==0).sum()),replace=True),rng.choice(np.flatnonzero(group==1),int((group==1).sum()),replace=True)]
   xx=np.column_stack([np.ones(len(ix)),group[ix],(score[ix]-score[ix].mean())/score[ix].std(ddof=1)])
   rec={'cohort':label,'bootstrap':b,'valid':False}
   if np.linalg.matrix_rank(xx)==3:
    try:zz=signed_z(Y[ix],xx)
    except ValueError:rec['invalid_reason']='unit_leverage_HC3_not_estimable';draws.append(rec);continue
    good=np.isfinite(zz)
    if good.all():
     single=stats.spearmanr(train['GSE63514'],zz).statistic;joint=stats.spearmanr(train['Joint'],zz).statistic
     if np.isfinite(single) and np.isfinite(joint):rec.update({'valid':True,'rho_single':single,'rho_joint':joint,'delta':joint-single});valid.append(joint-single)
   draws.append(rec)
  row.update({'bootstrap_attempts':1000,'valid_bootstraps':len(valid),'delta_ci_low':np.quantile(valid,.025),'delta_ci_high':np.quantile(valid,.975),'interval_interpretation':'unstable_small_n_descriptive_only' if len(cc)<10 or len(valid)<800 else 'exploratory_pointwise_conditional_interval'});summary.append(row);print(row,flush=True)
 pd.DataFrame(draws).to_csv(R/'discovery_strategy_bootstrap.tsv.gz',sep='\t',index=False,compression='gzip');pd.concat(full).to_csv(R/'new_study_genomewide_statistics.tsv.gz',sep='\t',index=False,compression='gzip');out=pd.DataFrame(summary);out.to_csv(R/'discovery_strategy_benchmark.tsv',sep='\t',index=False)
 fig,axs=plt.subplots(1,2,figsize=(11.5,3.9))
 for i,row in out.iterrows():
  axs[0].plot([0,1,2,3],[row.rho_GSE63514,row.rho_GSE138080,row.rho_GSE7803,row.rho_Joint],'-o',label=row.cohort)
  axs[1].errorbar(row.delta_joint_minus_single,i,xerr=[[row.delta_joint_minus_single-row.delta_ci_low],[row.delta_ci_high-row.delta_joint_minus_single]],fmt='o',capsize=3,color='#3a617b')
 axs[0].set(xticks=[0,1,2,3],xticklabels=['GSE63514','GSE138080','GSE7803','Joint'],ylabel='Rank correlation with new study',title='A  All individual and joint discovery rankings');axs[0].legend(frameon=False,fontsize=8)
 axs[1].axvline(0,color='grey',ls='--');axs[1].set(yticks=range(len(out)),yticklabels=out.cohort,xlabel='Joint minus single correlation\nPointwise 95% participant-bootstrap interval',title='B  Conditional uncertainty in transfer difference')
 fig.tight_layout();fig.savefig(F/'13_discovery_strategy_transfer.png',dpi=300);fig.savefig(F/'13_discovery_strategy_transfer.pdf');plt.close(fig)

if __name__=='__main__':main()

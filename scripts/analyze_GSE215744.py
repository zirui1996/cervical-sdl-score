import pandas as pd,numpy as np,statsmodels.api as sm,json
from statsmodels.stats.multitest import multipletests
from run_analysis import P,R
from extension_utils import extra_scores
counts=[];tpms=[]
for p in sorted(P.glob('GSM66381*_genes.deg_cx-swab*.tsv.gz')):
 d=pd.read_csv(p,sep='\t').dropna(subset=['gene']).groupby('gene')[['count','RPKM']].sum()
 counts.append(d['count'].rename(p.name[:10]));tpms.append(d.RPKM.rename(p.name[:10]))
ct=pd.concat(counts,axis=1).fillna(0);tp=pd.concat(tpms,axis=1).fillna(0)
assert ct.shape[1]==9 and (ct.values>=0).all() and np.allclose(ct,ct.round())
assert (tp.values>=0).all() and np.isfinite(tp.values).all()
rows=[]
for label,e in [('source_counts',np.log2(ct.div(ct.sum(),axis=1)*1e6+1)),('source_RPKM_same_patients',np.log2(tp+1))]:
 s,cov=extra_scores(e);s.to_csv(P/f'GSE215744_{label}_scores.tsv',sep='\t')
 ss=s.iloc[:6];g=np.array([0,0,0,0,1,1]);y=(ss.Joint_program-ss.Joint_program.mean())/ss.Joint_program.std()
 assert ss[['Joint_program','Composite']].notna().all().all()
 X=pd.DataFrame({'const':1.,'group':g,'Composite_z':(ss.Composite-ss.Composite.mean())/ss.Composite.std()},index=ss.index)
 f=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);ci=f.conf_int().loc['group']
 rows.append(dict(processing=label,n=6,HSIL_CIS=4,SQCC=2,beta=f.params['group'],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues['group'],residual_df=f.df_resid))
o=pd.DataFrame(rows);o['holm_two_processing_tests']=multipletests(o.p,method='holm')[1];o.to_csv(R/'GSE215744_fixed_program_transfer.tsv',sep='\t',index=False);print(o.to_string(index=False))

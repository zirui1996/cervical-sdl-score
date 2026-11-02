import io,gzip,tarfile,json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from run_analysis import RAW,P,R,design
from extension_utils import extra_scores

def main():
 text=(RAW/'GPL6699_metadata.txt').read_text().split('!platform_table_begin\n')[1].split('!platform_table_end')[0];a=pd.read_csv(io.StringIO(text),sep='\t',dtype=str)
 m=pd.read_csv(P/'GPL6699_probe_map_via_official_alternative.tsv',sep='\t',dtype=str).set_index('ID').gene
 a['position']=((a.MetaColumn.astype(int)-1)*4+a.MetaRow.astype(int)).astype(str)+'_'+a.Column+'_'+a.Row;ann=a.set_index('position')['ID'].map(m)
 ref=[]
 with tarfile.open(RAW/'GSE14404_RAW.tar') as t:
  for member in t.getmembers():
   if not member.name.startswith('GSM'):continue
   lines=gzip.decompress(t.extractfile(member).read()).decode('latin1').splitlines();i=next(i for i,l in enumerate(lines) if l.startswith('Block'))
   d=pd.read_csv(io.StringIO('\n'.join(lines[i:])),sep='\t',usecols=['Block','Column','Row','F532 Median'])
   ix=d.Block.astype(str)+'_'+d.Column.astype(str)+'_'+d.Row.astype(str);v=np.log2(d['F532 Median'].clip(lower=1));v.index=ix.map(ann);ref.append(v[v.index.notna()].groupby(level=0).median().rename(member.name.split('.')[0]))
 e=pd.concat(ref,axis=1);e.to_csv(P/'GSE14404_reference_Cy3_gene_log_intensity.tsv.gz',sep='\t',compression='gzip');rs,_=extra_scores(e);rs.to_csv(P/'GSE14404_reference_Cy3_scores.tsv',sep='\t')
 c=pd.read_csv(P/'GSE14404_independent_metadata.tsv',sep='\t',index_col=0);c=c.loc[c.primary_included].copy();c['batch_prefix']=c.description.str.lower().str.startswith('hoae').astype(int)
 assert pd.crosstab(c.group,c.batch_prefix).to_dict()=={0:{0:0,1:23},1:{0:3,1:4}}
 s=pd.read_csv(P/'GSE14404_raw_Cy5_scores.tsv',sep='\t',index_col=0)
 re=pd.read_csv(P/'GSE14404_gene_log_expression.tsv.gz',sep='\t',index_col=0);ratio,_=extra_scores(re)
 rows=[]
 def fit(name,cc,ss,ref_response=False,batch=False,ref_cov=False):
  X=pd.DataFrame({'const':1.,'group':cc.group.astype(float)},index=cc.index) if ref_response else design(cc,ss)
  if batch:X['batch_prefix']=cc.batch_prefix.astype(float)
  if ref_cov:X['reference_score_z']=(rs.loc[cc.index,'Joint_program']-rs.loc[cc.index,'Joint_program'].mean())/rs.loc[cc.index,'Joint_program'].std()
  v=(rs if ref_response else ss).loc[cc.index,'Joint_program'];y=(v-v.mean())/v.std();assert np.linalg.matrix_rank(X)==X.shape[1]
  f=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);ci=f.conf_int().loc['group'];rows.append(dict(analysis=name,n=len(cc),HSIL=int(cc.group.eq(0).sum()),Cancer=int(cc.group.eq(1).sum()),beta=f.params.group,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues.group,residual_df=f.df_resid,condition_number=float(np.linalg.cond(X))))
 fit('Cy5_plus_prefix_batch',c,s,batch=True)
 fit('Cy5_early_prefix_only',c.loc[c.batch_prefix.eq(1)],s)
 fit('Cy3_universal_reference_negative_control',c,s,ref_response=True)
 fit('Cy5_plus_reference_score',c,s,ref_cov=True)
 fit('log_ratio_plus_prefix_batch',c,ratio,batch=True)
 out=pd.DataFrame(rows);out['holm_5_assay_checks']=multipletests(out.p,method='holm')[1];out.to_csv(R/'GSE14404_assay_control_models.tsv',sep='\t',index=False);print(out.to_string(index=False))

if __name__=='__main__':main()

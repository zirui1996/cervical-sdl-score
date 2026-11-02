"""New independent two-channel tissue source: raw target intensity primary."""
import io,json,gzip,tarfile
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,RAW,P,R,F,expression,design
from extension_utils import read_meta,extra_scores
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 text=(RAW/'GPL6699_metadata.txt').read_text().split('!platform_table_begin\n')[1].split('!platform_table_end')[0]
 a=pd.read_csv(io.StringIO(text),sep='\t',dtype=str)
 m=pd.read_csv(P/'GPL6699_probe_map_via_official_alternative.tsv',sep='\t',dtype=str).set_index('ID').gene
 # This platform stores MetaColumn=1..12 and MetaRow=1..4; GPR blocks are column-major.
 a['block']=(a.MetaColumn.astype(int)-1)*4+a.MetaRow.astype(int)
 a['position']=a.block.astype(str)+'_'+a.Column+'_'+a.Row
 assert a.position.is_unique
 pos=a.set_index('position')['ID'];ann=a.set_index('position')['ID'].map(m)
 c=read_meta('GSE14404_series_matrix.txt.gz');c['stage']=np.select([c.title.str.startswith('CIN3'),c.title.str.startswith('CIN1'),c.title.str.startswith('NORMAL')],['HSIL','CIN1','Normal'],default='Cancer')
 c['group']=c.stage.eq('Cancer').astype(int);c['primary_included']=c.stage.eq('HSIL')|c.cancer_hpe.eq('SCC')
 assert c.primary_included.sum()==30 and c.loc[c.primary_included,'group'].value_counts().to_dict()=={1:27,0:3}
 c.to_csv(P/'GSE14404_independent_metadata.tsv',sep='\t')
 ratio=expression('GSE14404',m);raw=[];audits=[]
 with tarfile.open(RAW/'GSE14404_RAW.tar') as t:
  for member in t.getmembers():
   if not member.name.startswith('GSM'):continue
   sid=member.name.split('.')[0];assert sid in c.index
   text=gzip.decompress(t.extractfile(member).read()).decode('latin1');lines=text.splitlines();i=next(i for i,l in enumerate(lines) if 'F635 Median' in l and l.startswith('Block'))
   d=pd.read_csv(io.StringIO('\n'.join(lines[i:])),sep='\t',usecols=['Block','Column','Row','ID','F635 Median','F532 Median'])
   d['position']=d.Block.astype(str)+'_'+d.Column.astype(str)+'_'+d.Row.astype(str);assert d.position.is_unique
   d['native_id']=d.position.map(pos);d['gene']=d.position.map(ann)
   assert d.native_id.notna().all();assert d[['F635 Median','F532 Median']].ge(0).all().all()
   v=np.log2(d['F635 Median'].clip(lower=1));r=np.log2(d['F635 Median'].clip(lower=1)/d['F532 Median'].clip(lower=1))
   # Reconstruct author log-ratio ordering independently to check spatial probe/channel mapping.
   with gzip.open(RAW/'GSE14404_series_matrix.txt.gz','rt') as f:
    s=f.read().split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0]
   author=pd.read_csv(io.StringIO(s),sep='\t',index_col=0);author.index=author.index.astype(str)
   av=author[sid].reindex(d.native_id).to_numpy();valid=np.isfinite(av)&np.isfinite(r)
   rho=float(np.corrcoef(av[valid],r[valid])[0,1])
   assert rho>.95,(sid,rho)
   v.index=d.gene;v=v[v.index.notna()].groupby(level=0).median();raw.append(v.rename(sid))
   audits.append(dict(sample=sid,rows=len(d),gene_mapped=int(d.gene.notna().sum()),source_ratio_vs_reconstructed_correlation=rho,source_ratio_center_shift=float(np.nanmedian(av-r)),target_channel='F635 Median / Cy5'))
 e=pd.concat(raw,axis=1);assert set(e.columns)==set(c.index) and e.notna().all().all()
 e.to_csv(P/'GSE14404_raw_Cy5_gene_log_intensity.tsv.gz',sep='\t',compression='gzip')
 s,cov=extra_scores(e);s.to_csv(P/'GSE14404_raw_Cy5_scores.tsv',sep='\t');cov.to_csv(R/'GSE14404_score_coverage.tsv',sep='\t',index=False)
 rs,_=extra_scores(ratio);rows=[]
 for name,ss,mask in [('Cy5_primary_SCC',s,c.primary_included),('Cy5_include_unspecified_carcinoma',s,c.stage.isin(['HSIL','Cancer'])),('author_log_ratio_different_scale',rs,c.primary_included)]:
  cc=c.loc[mask];X=design(cc,ss);y=(ss.loc[cc.index,'Joint_program']-ss.loc[cc.index,'Joint_program'].mean())/ss.loc[cc.index,'Joint_program'].std()
  assert y.notna().all() and np.linalg.matrix_rank(X)==X.shape[1]
  f=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);ci=f.conf_int().loc['group']
  rows.append(dict(analysis=name,n=len(cc),HSIL=int(cc.group.eq(0).sum()),Cancer=int(cc.group.eq(1).sum()),beta=f.params.group,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues.group,residual_df=f.df_resid,condition_number=float(np.linalg.cond(X))))
 out=pd.DataFrame(rows);out['holm_3_program_tests']=multipletests(out.p,method='holm')[1];out.to_csv(R/'GSE14404_independent_fixed_program.tsv',sep='\t',index=False)
 pd.DataFrame(audits).to_csv(R/'GSE14404_raw_channel_probe_audit.tsv',sep='\t',index=False)
 (R/'GSE14404_independent_checks.json').write_text(json.dumps({'raw_samples':len(audits),'all_sample_and_physical_probe_identities_verified':True,'all_reconstructed_ratios_correlate_above_095':True,'target_channel':'Cy5 / F635 foreground median; no background subtraction','minimum_ratio_correlation':min(x['source_ratio_vs_reconstructed_correlation'] for x in audits),'unmapped_frozen_gene':'PINLYP','measured_frozen_genes':9,'no_feature_reselection':True,'not_RNAseq_validation':True},indent=2))
 fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
 for j,stage in enumerate(['Normal','CIN1','HSIL','Cancer']):
  vals=s.loc[c.index[c.stage.eq(stage)],'Joint_program'];ax.scatter(np.linspace(j-.12,j+.12,len(vals)),vals,s=25,color=plt.cm.viridis((j+.5)/4));ax.set(xticks=range(4),xticklabels=['Normal (5)','CIN1 (4)','CIN3 (3)','Carcinoma (28)'],ylabel='Frozen score | raw Cy5 abundance ranks',title='GSE14404 independent tissue source')
 fig.savefig(F/'22_GSE14404_independent_tissue.png',dpi=240);fig.savefig(F/'22_GSE14404_independent_tissue.pdf');plt.close(fig)
 print(out.to_string(index=False));print(cov.to_string(index=False))

if __name__=='__main__':main()

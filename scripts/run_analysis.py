"""Fixed proliferation-adjusted discovery; validation expression loaded after candidate lock."""
from pathlib import Path
import os,csv,gzip,hashlib,io,json,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];RAW=ROOT/'data/raw';P=ROOT/'data/processed';R=ROOT/'results';F=ROOT/'figures';SEED=20261029
for f in [P,R,F]:f.mkdir(parents=True,exist_ok=True)
FIXED=['MKI67','PCNA','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','TOP2A','CDK1','CCNB1','CCNB2']

def dump(name,data):
 (R/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,default=lambda x:x.item() if hasattr(x,'item') else str(x)))

def metadata(accession):
 m={};chars=[];series=[]
 with gzip.open(RAW/f'{accession}_series_matrix.txt.gz','rt') as handle:
  for line in handle:
   if line.startswith('!series_matrix_table_begin'):break
   if not line.startswith(('!Sample_','!Series_')):continue
   vals=next(csv.reader([line],delimiter='\t'));key=vals[0]
   if key.startswith('!Series_'):series.append({'field':key,'value':' | '.join(vals[1:])})
   elif key=='!Sample_characteristics_ch1':chars.append(vals[1:])
   elif key not in m:m[key]=vals[1:]
 c=pd.DataFrame(index=m['!Sample_geo_accession']);c.index.name='sample'
 for key in ['!Sample_title','!Sample_source_name_ch1','!Sample_data_processing','!Sample_extract_protocol_ch1']:
  if key in m:c[key.replace('!Sample_','')]=m[key]
 for i,vals in enumerate(chars):
  key=vals[0].split(':',1)[0] if ':' in vals[0] else f'characteristic_{i}'
  c[key]=[v.split(':',1)[1].strip() if ':' in v else v for v in vals]
 if accession=='GSE63514':
  c['stage']=c.title.str.split('-').str[0];c['whole_section']=c.dissection.eq('whole section').astype(int)
  c['primary_included']=c.stage.isin(['CIN3','Cancer'])&c.whole_section.eq(0)
  c['group']=c.stage.eq('Cancer').astype(int);c['excluded_test']=False
 elif accession=='GSE138080':
  c['stage']=np.select([c.title.str.startswith('Normal'),c.title.str.startswith('HGCIN'),c.title.str.startswith('CxSCC')],['Normal','HSIL','Cancer'],default='Unknown')
  c['group']=c.stage.eq('Cancer').astype(int);c['primary_included']=c.stage.isin(['HSIL','Cancer']);c['excluded_test']=False
 else:
  c['excluded_test']=c.title.str.startswith('Test Sample')
  c['stage']=np.select([c.title.str.startswith('Test Sample'),c.title.str.startswith('Normal'),c.title.str.startswith('high grade'),c.title.str.startswith('squamous cell')],['Test','Normal','HSIL','Cancer'],default='Unknown')
  c['group']=c.stage.eq('Cancer').astype(int);c['primary_included']=c.stage.isin(['HSIL','Cancer'])&~c.excluded_test
 c.to_csv(P/f'{accession}_metadata.tsv',sep='\t');pd.DataFrame(series).to_csv(P/f'{accession}_source_metadata.tsv',sep='\t',index=False)
 return c

def annotation(platform):
 with gzip.open(RAW/f'{platform}.annot.gz','rt') as h:
  lines=[line for line in h if line.strip() and not line.startswith(('!','#','^'))]
 a=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str)[['ID','Gene symbol']].dropna()
 a=a[a['Gene symbol'].str.fullmatch('[A-Za-z0-9_.-]+')&~a['Gene symbol'].eq('---')].drop_duplicates('ID')
 a.to_csv(P/f'{platform}_unambiguous_probe_map.tsv',sep='\t',index=False)
 return a.set_index('ID')['Gene symbol']

def expression(accession,ann):
 with gzip.open(RAW/f'{accession}_series_matrix.txt.gz','rt') as h:
  for line in h:
   if line.startswith('!series_matrix_table_begin'):break
  lines=[]
  for line in h:
   if line.startswith('!series_matrix_table_end'):break
   lines.append(line)
 e=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0);e.index=e.index.astype(str)
 rawq=np.nanquantile(e.values,[0,.01,.5,.99,1]);transform='none; author log-scale matrix'
 if rawq[3]>100 or rawq[4]-rawq[0]>50:
  assert rawq[0]>=0,'Negative values incompatible with automatic log transform'
  e=np.log2(e.clip(lower=1));transform='log2 of positive author-processed linear intensities; zeros clipped at 1'
 symbols=ann.reindex(e.index);probe_n=len(e);mapped_n=int(symbols.notna().sum())
 e=e.loc[symbols.notna()].copy();e.index=symbols.dropna().values;e=e.groupby(level=0,sort=True).median()
 e.to_csv(P/f'{accession}_gene_log_expression.tsv.gz',sep='\t',compression='gzip')
 dump(f'{accession}_processing.json',{'raw_value_quantiles':rawq.tolist(),'transformation':transform,'original_probes':probe_n,'unambiguously_mapped_probes':mapped_n,'genes':len(e),'samples':e.shape[1],'probe_collapse':'median without reference to stage'})
 return e

def gene_sets():
 sets={}
 for line in (RAW/'h.all.v2025.1.Hs.symbols.gmt').read_text().splitlines():
  v=line.split('\t');sets[v[0]]=set(v[2:])
 return {'E2F':sets['HALLMARK_E2F_TARGETS'],'G2M':sets['HALLMARK_G2M_CHECKPOINT'],'Fixed':set(FIXED),'Composite':sets['HALLMARK_E2F_TARGETS']|sets['HALLMARK_G2M_CHECKPOINT']}

def scores(e,common,sets):
 universe=sorted(set(common)&set(e.index));ranks=e.loc[universe].rank(axis=0,pct=True)
 out=pd.DataFrame(index=e.columns);coverage=[]
 for name,genes in sets.items():
  valid=sorted(set(universe)&genes)
  out[name]=ranks.loc[valid].mean(axis=0)
  coverage.append({'score':name,'defined_genes':len(genes),'available_genes':len(valid),'genes':';'.join(valid),'ranking_universe':len(universe)})
 return out,pd.DataFrame(coverage)

def design(c,s,kind='Composite',quadratic=False,dissection=False):
 v=s.loc[c.index,kind];z=(v-v.mean())/v.std(ddof=1)
 X=pd.DataFrame({'const':1.,'group':c.group.astype(float),'proliferation_z':z},index=c.index)
 if quadratic:X['proliferation_squared']=z**2
 if dissection:X['whole_section']=c.whole_section.astype(float)
 return X

def ols_all(e,X):
 """Vectorized OLS HC3 with t residual-df inference; verified against statsmodels."""
 Y=e[X.index].T.values.astype(float);A=X.values.astype(float);inv=np.linalg.inv(A.T@A);B=inv@A.T@Y
 residual=Y-A@B;leverage=np.einsum('ij,jk,ik->i',A,inv,A)
 w=(A@inv)[:,1];variance=np.sum((residual/(1-leverage[:,None]))**2*w[:,None]**2,axis=0)
 se=np.sqrt(variance);t=B[1]/se;df=len(X)-X.shape[1];p=2*stats.t.sf(abs(t),df);crit=stats.t.ppf(.975,df)
 d=pd.DataFrame({'gene':e.index,'beta':B[1],'se_HC3':se,'t':t,'p':p,'ci_low':B[1]-crit*se,'ci_high':B[1]+crit*se,'residual_df':df})
 d['fdr']=multipletests(np.nan_to_num(d.p,nan=1),method='fdr_bh')[1]
 return d.sort_values(['p','gene']).reset_index(drop=True)

def auc_bootstrap(values,y,direction):
 rng=np.random.default_rng(SEED);a=direction*values[y==1];b=direction*values[y==0]
 aa=rng.choice(a,(2000,len(a)),replace=True);bb=rng.choice(b,(2000,len(b)),replace=True)
 rr=stats.rankdata(np.column_stack([aa,bb]),axis=1);aucs=(rr[:,:len(a)].sum(axis=1)-len(a)*(len(a)+1)/2)/(len(a)*len(b))
 r=stats.rankdata(np.r_[a,b]);auc=(r[:len(a)].sum()-len(a)*(len(a)+1)/2)/(len(a)*len(b))
 return {'unadjusted_oriented_auc':auc,'auc_ci_low':np.quantile(aucs,.025),'auc_ci_high':np.quantile(aucs,.975)}

def diagnostics(label,c,s):
 X=design(c,s);rho=stats.pointbiserialr(c.group,s.loc[c.index,'Composite']).statistic
 groups=[s.loc[c.index[c.group.eq(k)],'Composite'] for k in [0,1]];low=max(g.min() for g in groups);high=min(g.max() for g in groups)
 return {'cohort':label,'n':len(c),'precancer':int(c.group.eq(0).sum()),'cancer':int(c.group.eq(1).sum()),'proliferation_group_r':rho,'group_proliferation_VIF':1/(1-rho*rho),'design_condition_number':np.linalg.cond(X.values),'common_support_min':low,'common_support_max':high,'precancer_in_common_support':int(groups[0].between(low,high).sum()),'cancer_in_common_support':int(groups[1].between(low,high).sum())}

def sensitivities(label,c,e,s,candidates):
 rows=[]
 def add(name,cc,XX):
  if len(cc)<6 or cc.group.value_counts().min()<3:return
  dd=ols_all(e.loc[candidates],XX)
  for row in dd.to_dict('records'):row.update({'cohort':label,'sensitivity':name,'n':len(cc),'cancer':int(cc.group.sum())});rows.append(row)
 for kind in ['E2F','G2M','Fixed']:add(kind,c,design(c,s,kind))
 add('quadratic_proliferation',c,design(c,s,quadratic=True))
 lo=max(s.loc[c.index[c.group.eq(k)],'Composite'].min() for k in [0,1]);hi=min(s.loc[c.index[c.group.eq(k)],'Composite'].max() for k in [0,1]);sub=c.loc[s.loc[c.index,'Composite'].between(lo,hi)]
 if len(sub.group.value_counts())==2 and sub.group.value_counts().min()>=5:add('common_support_only',sub,design(sub,s))
 for sample in c.index:
  sub=c.drop(index=sample);add('leave_out_'+sample,sub,design(sub,s))
 return rows

def main():
 platforms={'GSE63514':'GPL570','GSE138080':'GPL4133','GSE7803':'GPL96'}
 meta={k:metadata(k) for k in platforms};ann={k:annotation(p) for k,p in platforms.items()}
 common=set.intersection(*(set(a.values) for a in ann.values()));sets=gene_sets();exclude=set.union(*(v for v in sets.values()))
 e=expression('GSE63514',ann['GSE63514']);sc,cov=scores(e,common,sets);sc.to_csv(P/'GSE63514_proliferation_scores.tsv',sep='\t');cov.to_csv(R/'GSE63514_score_coverage.tsv',sep='\t',index=False)
 c=meta['GSE63514'].loc[meta['GSE63514'].primary_included].copy()
 pool=e[c.index];eligible=pool.index.isin(common-exclude)&pool.notna().all(axis=1)&pool.std(axis=1).gt(0)&pool.median(axis=1).ge(4)&(pool.quantile(.75,axis=1)-pool.quantile(.25,axis=1)).ge(.25)
 test=pool.loc[eligible];test.index.name='gene';pd.Series(sorted(common),name='gene').to_csv(R/'platform_common_gene_universe.tsv',sep='\t',index=False)
 pd.Series(test.index,name='gene').to_csv(R/'primary_tested_gene_universe.tsv',sep='\t',index=False)
 adj=ols_all(test,design(c,sc));unadj=ols_all(test,pd.DataFrame({'const':1.,'group':c.group.astype(float)},index=c.index))
 adj=adj.merge(unadj[['gene','beta','p','fdr']],on='gene',suffixes=('','_unadjusted'))
 adj.to_csv(R/'discovery_all_genes.tsv',sep='\t',index=False)
 qualified=adj[adj.fdr.lt(.05)&adj.beta.abs().ge(.5)];selected=qualified.head(5).copy();genes=selected.gene.tolist()
 lock={'locked_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'discovery_n':len(c),'discovery_CIN3':int(c.group.eq(0).sum()),'discovery_cancer_laser_capture':int(c.group.eq(1).sum()),'tested_genes':len(test),'adjusted_FDR_significant':int(adj.fdr.lt(.05).sum()),'qualified_candidates':len(qualified),'candidates':selected.to_dict('records'),'plan_sha256':hashlib.sha256((ROOT/'ANALYSIS_PLAN_v1.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'score_definitions':{k:sorted(v) for k,v in sets.items()}}
 lockfile=R/'candidate_lock.json'
 if lockfile.exists():
  prior=json.loads(lockfile.read_text());assert [r['gene'] for r in prior['candidates']]==genes,'Candidate lock mismatch; do not overwrite.'
 else:dump('candidate_lock.json',lock)
 print('LOCK',json.dumps({k:v for k,v in lock.items() if k!='score_definitions'}),flush=True)
 validation=[];sens=[];diags=[diagnostics('GSE63514',c,sc)];cohorts={'GSE63514':(c,e,sc)}
 if genes:
  sens+=sensitivities('GSE63514',c,e,sc,genes)
  # Discovery-wide section and precursor-grade sensitivities retain the frozen target universe and candidate lock.
  for name,mask in [('all_sections_with_dissection_adjustment',meta['GSE63514'].stage.isin(['CIN3','Cancer'])),('CIN2_CIN3_laser_only',meta['GSE63514'].stage.isin(['CIN2','CIN3','Cancer'])&meta['GSE63514'].whole_section.eq(0))]:
   sub=meta['GSE63514'].loc[mask];XX=design(sub,sc,dissection=(name=='all_sections_with_dissection_adjustment'))
   res=ols_all(e.loc[genes],XX)
   for row in res.to_dict('records'):row.update({'cohort':'GSE63514','sensitivity':name,'n':len(sub),'cancer':int(sub.group.sum())});sens.append(row)
 # External expression is opened for the first time only after the candidate lock above.
 for label in ['GSE138080','GSE7803']:
  ve=expression(label,ann[label]);vs,vcov=scores(ve,common,sets);vs.to_csv(P/f'{label}_proliferation_scores.tsv',sep='\t');vcov.to_csv(R/f'{label}_score_coverage.tsv',sep='\t',index=False)
  cc=meta[label].loc[meta[label].primary_included].copy();diags.append(diagnostics(label,cc,vs));cohorts[label]=(cc,ve,vs)
  for row in selected.to_dict('records'):
   gene=row['gene'];rr={'cohort':label,'gene':gene,'discovery_beta':row['beta'],'discovery_fdr':row['fdr']}
   if gene not in ve.index:rr.update({'p':1.,'status':'unmeasured','direction_concordant':False});validation.append(rr);continue
   ccg=cc.loc[ve.loc[gene,cc.index].notna()];fit=ols_all(ve.loc[[gene]],design(ccg,vs)).iloc[0].to_dict();fit.pop('fdr');rr.update(fit);rr.update({'n':len(ccg),'precancer':int(ccg.group.eq(0).sum()),'cancer':int(ccg.group.eq(1).sum()),'direction_concordant':bool(np.sign(rr['beta'])==np.sign(row['beta'])),'status':'tested'})
   rr.update(auc_bootstrap(ve.loc[gene,ccg.index].values,ccg.group.values,np.sign(row['beta'])));validation.append(rr)
  if genes:sens+=sensitivities(label,cc,ve,vs,[g for g in genes if g in ve.index])
 val=pd.DataFrame(validation)
 if len(val):
  val['holm_p']=val.groupby('cohort').p.transform(lambda x:multipletests(x,method='holm')[1]);val['replicated']=val.holm_p.lt(.05)&val.direction_concordant
  val.to_csv(R/'external_validation.tsv',sep='\t',index=False)
 pd.DataFrame(sens).to_csv(R/'sensitivity_analyses.tsv',sep='\t',index=False);pd.DataFrame(diags).to_csv(R/'design_diagnostics.tsv',sep='\t',index=False)
 # Numerical validation of the vectorized inferential engine against an independent library implementation.
 checks=[]
 for gene in [adj.gene.iloc[0],adj.gene.iloc[len(adj)//2],adj.gene.iloc[-1]]:
  fit=sm.OLS(test.loc[gene],design(c,sc)).fit(cov_type='HC3',use_t=True);rr=adj.set_index('gene').loc[gene]
  assert np.isclose(fit.params['group'],rr.beta) and np.isclose(fit.pvalues['group'],rr.p,rtol=1e-7,atol=1e-12);checks.append(gene)
 dump('numerical_verification.json',{'HC3_matches_statsmodels_use_t':checks,'passed':True})
 both=[]
 if len(val):both=[g for g in genes if val.loc[val.gene.eq(g),'replicated'].sum()==2]
 summary={'primary_discovery':{'CIN3':40,'cancer_laser':20,'whole_section_excluded':8},'primary_validation':{'HSIL':15,'SCC':10},'secondary_validation':{'HSIL':7,'SCC':21,'technical_tests_excluded':3},'tested_genes':len(test),'discovery_FDR_significant':int(adj.fdr.lt(.05).sum()),'qualified_candidates':len(qualified),'locked_genes':genes,'replicated_primary':val.loc[val.cohort.eq('GSE138080')&val.replicated,'gene'].tolist() if len(val) else [],'replicated_secondary':val.loc[val.cohort.eq('GSE7803')&val.replicated,'gene'].tolist() if len(val) else [],'replicated_both':both,'caution':'Conditional cross-sectional association; statistical adjustment does not prove biological independence, causality or future progression.'}
 dump('analysis_summary.json',summary);figures(adj,cohorts,genes,val);print('SUMMARY',json.dumps(summary),flush=True)
 if len(val):print(val.to_string(index=False),flush=True)

def figures(adj,cohorts,genes,val):
 plt.rcParams.update({'font.size':10,'figure.dpi':150,'axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(1,3,figsize=(14,4))
 for ax,(label,(c,e,s)) in zip(axs,cohorts.items()):
  values=[s.loc[c.index[c.group.eq(k)],'Composite'] for k in [0,1]];ax.boxplot(values,tick_labels=['CIN3/HSIL','Invasive SCC'],showfliers=False)
  rng=np.random.default_rng(SEED)
  for i,vs in enumerate(values):ax.scatter(rng.normal(i+1,.04,len(vs)),vs,s=15,alpha=.6)
  ax.set(title=label,ylabel='Mean fractional-rank proliferation score')
 fig.suptitle('Proliferation score: prespecified E2F + G2M markers');fig.tight_layout();fig.savefig(F/'01_proliferation_distributions.png');fig.savefig(F/'01_proliferation_distributions.pdf');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,5));ax.scatter(adj.beta,-np.log10(adj.p.clip(1e-300)),s=7,c=np.where(adj.fdr<.05,'#bd5a44','#b6c3ce'),alpha=.6)
 for gene in genes:
  r=adj.set_index('gene').loc[gene];ax.annotate(gene,(r.beta,-np.log10(r.p)),fontsize=8)
 ax.set(xlabel='Invasive-group coefficient, proliferation adjusted',ylabel='-log10(HC3 P)',title='Laser-capture discovery: all tested genes');fig.tight_layout();fig.savefig(F/'02_adjusted_discovery.png');fig.savefig(F/'02_adjusted_discovery.pdf');plt.close(fig)
 if not genes:return
 fig,axs=plt.subplots(len(genes),3,figsize=(13,2.6*len(genes)),squeeze=False)
 for i,gene in enumerate(genes):
  for ax,(label,(c,e,s)) in zip(axs[i],cohorts.items()):
   if gene not in e.index:continue
   for k,color in [(0,'#457b9d'),(1,'#d56c4a')]:
    ix=c.index[c.group.eq(k)];ax.scatter(s.loc[ix,'Composite'],e.loc[gene,ix],s=18,color=color,label='SCC' if k else 'CIN3/HSIL',alpha=.8)
   ax.set(title=f'{gene} | {label}',xlabel='Proliferation rank score',ylabel='Author log2 expression')
   if i==0:ax.legend(fontsize=7)
 fig.tight_layout();fig.savefig(F/'03_locked_candidate_scatter.png');fig.savefig(F/'03_locked_candidate_scatter.pdf');plt.close(fig)
 fig,ax=plt.subplots(figsize=(9,7));rows=[]
 for gene in genes:
  rr=adj.set_index('gene').loc[gene];rows.append((gene+' / discovery',rr.beta,rr.ci_low,rr.ci_high))
  for _,rr in val[val.gene.eq(gene)].iterrows():
   if rr.status=='tested':rows.append((gene+' / '+rr.cohort,rr.beta,rr.ci_low,rr.ci_high))
 for i,(label,b,lo,hi) in enumerate(rows):ax.errorbar(b,i,xerr=np.array([[b-lo],[hi-b]]),fmt='o',capsize=3,color='#38566e')
 ax.set(yticks=range(len(rows)),yticklabels=[r[0] for r in rows],xlabel='Adjusted invasive-group log2 expression coefficient (95% CI)',title='Frozen candidates across independent cohorts');ax.axvline(0,color='grey',ls='--');ax.invert_yaxis();fig.tight_layout();fig.savefig(F/'04_cross_cohort_effects.png');fig.savefig(F/'04_cross_cohort_effects.pdf');plt.close(fig)

if __name__=='__main__':main()

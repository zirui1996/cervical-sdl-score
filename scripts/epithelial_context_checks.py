import gc,gzip,json
import numpy as np
import pandas as pd
from scipy.io import mmread
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from run_analysis import P,R,RAW,design,ols_all
from extension_utils import common_universe,extra_scores,frozen_targets
from analyze_stage_and_composition import load_all,responses

SETS={'Glandular_secretory':['PAX8','TFF3','AGR2','MUC5B','MUC5AC','KRT7'],'Squamous_basal':['KRT5','KRT14','TP63','DSG3','DSC3','SFN'],'Squamous_differentiated':['KRT4','KRT13','IVL','SPRR1A','SPRR3','KRT1','DSG1']}

def main():
 cohorts=load_all();old,meta,directions=frozen_targets();targets=list(directions);rows=[];diag=[];coverage=[]
 for label in ['GSE63514','GSE138080','GSE7803','GSE166466','GSE147009','GSE149763','GSE223804']:
  c,e,s=cohorts[label];cc=c.loc[c.primary_included];rank=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True);Y=responses(cc,e,s,targets)
  for name,genes in SETS.items():
   available=rank.index.intersection(genes);v=rank.loc[available].mean(axis=0);coverage.append({'cohort':label,'proxy':name,'defined':len(genes),'measured':len(available),'genes':';'.join(available)})
   if len(available)<3 or v.std()==0:continue
   X=design(cc,s);X[name+'_z']=(v-v.mean())/v.std(ddof=1);dd=ols_all(Y,X);dd['cohort']=label;dd['proxy']=name;dd['n']=len(cc);dd['interpretation']='post_hoc_context_sensitivity';rows.append(dd)
   diag.append({'cohort':label,'proxy':name,'n':len(cc),'residual_df':len(cc)-X.shape[1],'group_VIF':1/(1-sm.OLS(X.group,X.drop(columns='group')).fit().rsquared),'condition_number':np.linalg.cond(X)})
 out=pd.concat(rows,ignore_index=True);out['epithelial_context_fdr']=multipletests(out.p,method='fdr_bh')[1];out.to_csv(R/'epithelial_context_proxy_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(diag).to_csv(R/'epithelial_context_design_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(R/'epithelial_context_proxy_coverage.tsv',sep='\t',index=False)
 # Each sample is loaded in isolation; subtype annotation uses no target genes.
 cm=pd.read_csv(P/'GSE208653_cell_QC_marker_compartments.tsv.gz',sep='\t');pb=[];meta_rows=[];marker_summary=[]
 for f in sorted(RAW.glob('GSM636068*.matrix.mtx.gz')):
  donor=f.name.split('_')[0];base=f.name[:-len('.matrix.mtx.gz')];cells=cm[(cm.donor==donor)&cm.marker_compartment.eq('Epithelial')];stage=cells.stage.iloc[0]
  features=pd.read_csv(RAW/(base+'.features.tsv.gz'),sep='\t',header=None);symbols=features.iloc[:,1].astype(str).to_numpy();bc=pd.read_csv(RAW/(base+'.barcodes.tsv.gz'),sep='\t',header=None).iloc[:,0].astype(str);ix=pd.Index(bc).get_indexer(cells.barcode);assert (ix>=0).all()
  with gzip.open(f,'rb') as h:X=mmread(h).tocsr().astype(np.int32)
  X=X[:,ix].tocsr();lib=np.asarray(X.sum(axis=0)).ravel();scores=[];detection=[]
  for name,markers in SETS.items():
   vv=[];dd=[]
   for gene in markers:
    values=np.asarray(X[symbols==gene].sum(axis=0)).ravel();vv.append(np.log1p(values/lib*1e4));dd.append(values>0)
   scores.append(np.mean(vv,axis=0));detection.append(np.sum(dd,axis=0))
  score=np.array(scores);det=np.array(detection);order=np.argsort(score,axis=0);best=order[-1];second=order[-2];k=np.arange(len(ix));good=(score[best,k]>=.3)&(det[best,k]>=2)&((score[best,k]-score[second,k])>=.15);names=np.array(list(SETS),dtype=object)[best];names[~good]='Mixed_unspecified'
  for state in sorted(set(names)):
   mask=names==state;counts=np.asarray(X[:,mask].sum(axis=1)).ravel();vec=pd.DataFrame({'gene':symbols,'count':counts}).groupby('gene')['count'].sum();key=donor+'|'+state;pb.append(vec.rename(key));meta_rows.append({'key':key,'donor':donor,'stage':stage,'marker_state':state,'cells':int(mask.sum()),'fraction_of_epithelial':float(mask.mean())})
   for j,name in enumerate(SETS):marker_summary.append({'donor':donor,'stage':stage,'state':state,'marker_set':name,'mean_marker_score':float(score[j,mask].mean())})
  print(donor,pd.Series(names).value_counts().to_dict(),flush=True);del X;gc.collect()
 counts=pd.concat(pb,axis=1).fillna(0).astype(np.int64);counts.to_csv(P/'GSE208653_epithelial_state_pseudobulk_counts.tsv.gz',sep='\t',compression='gzip');e=np.log2(counts.div(counts.sum(),axis=1)*1e6+1);s,cov=extra_scores(e);m=pd.DataFrame(meta_rows).set_index('key').join(s[['Joint_program','Composite']]);m.loc[m.cells<30,'Joint_program']=np.nan;m.reset_index().to_csv(R/'GSE208653_epithelial_marker_states.tsv',sep='\t',index=False);pd.DataFrame(marker_summary).to_csv(R/'GSE208653_epithelial_state_marker_review.tsv',sep='\t',index=False)
 print('\nNew array program with each epithelial-context proxy:\n',out[out.cohort.eq('GSE166466')&out.gene.eq('Joint_program')].to_string(index=False));print('\nEpithelial state pseudobulk:\n',m[m.stage.isin(['HSIL','SCC'])].to_string())

if __name__=='__main__':main()

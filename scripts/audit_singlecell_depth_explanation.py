"""Post hoc fixed-score technical explanation check, no biological replication claims."""
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from run_analysis import ROOT,P,R,F
from extension_utils import common_universe
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 rng=np.random.default_rng(20261029);common=common_universe()
 lock=json.loads((R/'joint_program_lock.json').read_text());genes=[x['gene'] for x in lock['candidates']];directions=np.array([x['direction'] for x in lock['candidates']])
 summaries=[];replicates=[];contributions=[];samples=[]
 for scope in ['Epithelial','Squamous_basal']:
  objects=[]
  for study in ['GSE208653','EMTAB12305']:
   state=scope!='Epithelial';name='epithelial_marker_states' if state else 'specimen_compartment_programs'
   meta=pd.read_csv(R/f'{study}_{name}.tsv',sep='\t')
   if 'pseudobulk_key' in meta:meta=meta.rename(columns={'pseudobulk_key':'key'})
   if 'donor' in meta:meta=meta.rename(columns={'donor':'specimen'})
   part=meta.loc[meta.stage.isin(['HSIL','SCC'])&meta.cells.ge(30)&meta['marker_state' if state else 'compartment'].eq(scope)].copy()
   ct=pd.read_csv(P/f'{study}_{"epithelial_state" if state else "compartment"}_pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
   universe=[x for x in common if x in ct.index];assert set(genes)<=set(universe)
   loc=[universe.index(g) for g in genes];ct=ct.loc[:,part.key];assert (ct.values>=0).all() and np.allclose(ct,ct.round())
   for _,row in part.iterrows():
    vals=ct.loc[universe,row.key].to_numpy(dtype=np.int64);total=int(ct[row.key].sum());counts=np.r_[vals,total-int(vals.sum())]
    original=rankdata(vals)/len(vals);score=float((original[loc]*directions).mean())
    assert np.isclose(score,row.Joint_program,atol=1e-10),(study,scope,row.key,score,row.Joint_program)
    objects.append(dict(study=study,scope=scope,specimen=row.specimen,stage=row.stage,cells=int(row.cells),total=total,
                        counts=counts,loc=loc,universe=len(vals),original_score=score,original_components=original[loc]*directions))
  depth=min(x['total'] for x in objects);assert len(objects)==10
  sims=[]
  for obj in objects:
   draws=rng.multivariate_hypergeometric(obj['counts'],depth,size=200)
   assert (draws.sum(axis=1)==depth).all() and (draws<=obj['counts']).all()
   ranks=rankdata(draws[:,:-1],axis=1)/obj['universe'];scores=(ranks[:,obj['loc']]*directions).mean(axis=1)
   sims.append(scores)
   samples.append({k:obj[k] for k in ['study','scope','specimen','stage','cells','total','universe','original_score']}|dict(matched_depth=depth,rarefied_score_mean=float(scores.mean()),rarefied_score_sd=float(scores.std()),original_detected_common=int((obj['counts'][:-1]>0).sum()),rarefied_detected_common_mean=float((draws[:,:-1]>0).sum(axis=1).mean())))
  sims=np.array(sims)
  for study in ['GSE208653','EMTAB12305']:
   hi=np.array([x['study']==study and x['stage']=='HSIL' for x in objects]);ca=np.array([x['study']==study and x['stage']=='SCC' for x in objects]);base=np.array([x['original_score'] for x in objects]);delta=sims[ca].mean(axis=0)-sims[hi].mean(axis=0)
   original=float(base[ca].mean()-base[hi].mean())
   summaries.append(dict(scope=scope,study=study,HSIL=int(hi.sum()),SCC=int(ca.sum()),matched_UMIs=depth,original_SCC_minus_HSIL=original,rarefied_delta_mean=float(delta.mean()),technical_q025=float(np.quantile(delta,.025)),technical_q975=float(np.quantile(delta,.975)),fraction_technical_draws_positive=float((delta>0).mean()),not_patient_confidence_interval=True))
   for b,v in enumerate(delta):replicates.append(dict(scope=scope,study=study,replicate=b,delta=v))
   comp=np.array([x['original_components'] for x in objects]);diff=comp[ca].mean(axis=0)-comp[hi].mean(axis=0);assert np.isclose(diff.mean(),original)
   for g,v in zip(genes,diff):contributions.append(dict(scope=scope,study=study,gene=g,signed_rank_contribution_to_score_difference=v/len(genes)))
 out=pd.DataFrame(summaries);out.to_csv(R/'singlecell_depth_matched_summary.tsv',sep='\t',index=False)
 pd.DataFrame(samples).to_csv(R/'singlecell_depth_sample_diagnostics.tsv',sep='\t',index=False)
 pd.DataFrame(replicates).to_csv(R/'singlecell_depth_technical_replicates.tsv',sep='\t',index=False)
 pd.DataFrame(contributions).to_csv(R/'singlecell_direction_gene_contributions.tsv',sep='\t',index=False)
 (R/'mainline_depth_audit_checks.json').write_text(json.dumps({'seed':20261029,'replicates':200,'specimen_scopes':len(samples),'original_scores_exactly_reproduced':True,'all_draw_totals_and_without_replacement_bounds_verified':True,'gene_contributions_sum_to_original_difference':True,'plan_sha256':hashlib.sha256((ROOT/'MAINLINE_DECISION_AUDIT_v1.md').read_bytes()).hexdigest(),'not_independent_validation':True},indent=2))
 fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 for ax,scope in zip(axs,['Epithelial','Squamous_basal']):
  part=out[out.scope.eq(scope)]
  for i,(_,r) in enumerate(part.iterrows()):
   ax.scatter(i-.12,r.original_SCC_minus_HSIL,color='#ac684a',label='Original' if i==0 else None)
   ax.errorbar(i+.12,r.rarefied_delta_mean,yerr=[[r.rarefied_delta_mean-r.technical_q025],[r.technical_q975-r.rarefied_delta_mean]],fmt='o',capsize=4,color='#436c83',label='Depth matched (technical 95% range)' if i==0 else None)
  ax.axhline(0,color='grey',ls='--');ax.set(xticks=[0,1],xticklabels=['GSE208653','E-MTAB-12305'],ylabel='Mean fixed score: SCC minus HSIL',title=scope)
 axs[0].legend(fontsize=7);fig.savefig(F/'21_singlecell_depth_explanation.png',dpi=240);fig.savefig(F/'21_singlecell_depth_explanation.pdf');plt.close(fig)
 print(out.to_string(index=False))

if __name__=='__main__':main()

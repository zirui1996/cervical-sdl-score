import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,P,R,F,design,ols_all,dump
from extension_utils import extra_scores,frozen_targets,fit_targets,program_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_all():
 cohorts={}
 for label in ['GSE63514','GSE138080','GSE7803']:
  c=pd.read_csv(P/f'{label}_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/f'{label}_gene_log_expression.tsv.gz',sep='\t',index_col=0);s,cov=extra_scores(e);cohorts[label]=(c,e,s)
 for label in ['GSE166466','GSE149763','GSE147009','GSE51993','GSE75132']:
  c=pd.read_csv(P/f'{label}_extension_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/f'{label}_extension_gene_log_expression.tsv.gz',sep='\t',index_col=0);s=pd.read_csv(P/f'{label}_extension_scores.tsv',sep='\t',index_col=0);cohorts[label]=(c,e,s)
 c=pd.read_csv(P/'GSE223804_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/'GSE223804_log2_CPMplus1.tsv.gz',sep='\t',index_col=0);s,cov=extra_scores(e);cohorts['GSE223804']=(c,e,s)
 return cohorts

def responses(c,e,s,targets):
 Y=e.loc[e.index.intersection(targets),c.index].copy();v=s.loc[c.index,'Joint_program']
 if v.notna().all() and v.std()>0:Y.loc['Joint_program']=(v-v.mean())/v.std(ddof=1)
 return Y

def contrast_rows(label,c,e,s,targets,levels):
 cc=c.loc[c.stage.isin(levels)].copy();stage=cc.stage.copy()
 if label=='GSE63514':cc=cc.loc[cc.whole_section.eq(0)];stage=cc.stage
 X=pd.DataFrame({'const':1.,'level_1':stage.eq(levels[1]).astype(float),'level_2':stage.eq(levels[2]).astype(float)},index=cc.index)
 v=s.loc[cc.index,'Composite'];X['proliferation_z']=(v-v.mean())/v.std(ddof=1)
 contrasts=[(levels[1]+'_vs_'+levels[0],[0,1,0,0]),(levels[2]+'_vs_'+levels[1],[0,-1,1,0]),('late_minus_early',[0,-2,1,0])]
 Y=responses(cc,e,s,targets);rows=[]
 for gene in targets+['Joint_program']:
  for name,vector in contrasts:
   rr={'cohort':label,'target':gene,'contrast':name,'baseline':levels[0],'intermediate':levels[1],'late':levels[2],'n':len(cc)}
   if gene not in Y.index or Y.loc[gene].isna().any() or Y.loc[gene].std()==0:rr.update({'p':1.,'status':'nonestimable'});rows.append(rr);continue
   fit=sm.OLS(Y.loc[gene],X).fit(cov_type='HC3',use_t=True);ct=fit.t_test(vector);ci=np.asarray(ct.conf_int()).ravel()
   rr.update({'beta':float(np.asarray(ct.effect).ravel()[0]),'p':float(np.asarray(ct.pvalue)),'ci_low':float(ci[0]),'ci_high':float(ci[1]),'status':'tested'});rows.append(rr)
 return rows

def main():
 cohorts=load_all();old,metagenes,directions=frozen_targets();targets=list(directions);stage=[];composition=[];diag=[]
 for label in ['GSE63514','GSE138080','GSE7803','GSE166466','GSE147009','GSE51993']:
  c,e,s=cohorts[label]
  levels=['Normal','CIN3' if label=='GSE63514' else 'HSIL','Cancer']
  if label=='GSE147009':levels=['HPV16_normal','HSIL','Cancer']
  if label=='GSE51993':levels=['Normal','CIN1','HSIL']
  stage+=contrast_rows(label,c,e,s,targets,levels)
 staged=pd.DataFrame(stage);staged['extension_stage_fdr']=multipletests(staged.p,method='fdr_bh')[1];staged.to_csv(R/'stage_context_all_contrasts.tsv',sep='\t',index=False)
 # Same primary samples, add proxies without claiming calibrated cell fractions.
 for label in ['GSE63514','GSE138080','GSE7803','GSE166466','GSE147009','GSE149763','GSE223804']:
  c,e,s=cohorts[label];cc=c.loc[c.primary_included];XX=design(cc,s)
  for proxy in ['Immune_proxy','Stroma_proxy']:
   v=s.loc[cc.index,proxy];assert v.notna().all() and v.std()>0
   XX[proxy+'_z']=(v-v.mean())/v.std(ddof=1)
  Y=responses(cc,e,s,targets);rr=ols_all(Y,XX);rr['cohort']=label;rr['model']='Composite_plus_immune_stroma_proxies';rr['n']=len(cc);composition.append(rr)
  vif=1/(1-sm.OLS(XX.group,XX.drop(columns='group')).fit().rsquared)
  diag.append({'cohort':label,'n':len(cc),'residual_df':len(cc)-XX.shape[1],'condition_number':float(np.linalg.cond(XX)),'group_VIF':float(vif),'limits':'expression proxies, not fractions or causal controls'})
 composition=pd.concat(composition,ignore_index=True);composition['composition_extension_fdr']=multipletests(composition.p,method='fdr_bh')[1]
 composition['interpretation']=np.where(composition.residual_df.le(1),'not_interpretable_df1_with_high_collinearity','supportive_expression_proxy_model')
 composition.to_csv(R/'composition_proxy_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(diag).to_csv(R/'composition_design_diagnostics.tsv',sep='\t',index=False)
 # For the longitudinal endpoint, correct the union and program together.
 c,e,s=cohorts['GSE75132'];clinical=pd.DataFrame(fit_targets('GSE75132',c,e,s,targets,directions,True)+[program_fit('GSE75132',c,s,True)]);clinical['holm_union_and_program']=multipletests(clinical.p,method='holm')[1];clinical.to_csv(R/'prospective_endpoint_transfer.tsv',sep='\t',index=False)
 # Point summaries and development/new/RNA roles stay visibly separate.
 primary=[]
 for label in ['GSE63514','GSE138080','GSE7803','GSE166466','GSE149763','GSE147009','GSE223804']:
  c,e,s=cohorts[label];r=program_fit(label,c,s,True);r['role']='joint_discovery' if label in ['GSE63514','GSE138080','GSE7803'] else 'post_hoc_existing_RNA' if label=='GSE223804' else 'new_held_out_tissue';primary.append(r)
 primary=pd.DataFrame(primary);primary.to_csv(R/'joint_program_all_tissue_cohorts.tsv',sep='\t',index=False)
 fig,ax=plt.subplots(figsize=(8.5,5.1))
 for i,row in primary.iterrows():
  color='#627b8e' if row.role=='joint_discovery' else '#ad704a' if row.role=='new_held_out_tissue' else '#8d8d8d'
  ax.errorbar(row.beta,i,xerr=[[row.beta-row.ci_low],[row.ci_high-row.beta]],fmt='o',capsize=3,color=color)
 display_roles=primary.role.map({'joint_discovery':'Development','new_held_out_tissue':'New tissue','post_hoc_existing_RNA':'Prior RNA context'})
 ax.axvline(0,color='grey',ls='--');ax.set(yticks=range(len(primary)),yticklabels=primary.cohort+' | '+display_roles,xlabel='Adjusted program coefficient (SD units; 95% HC3 CI)',title='Fixed program across tissue cohorts');ax.invert_yaxis();fig.tight_layout();fig.savefig(F/'10_program_cross_cohort.png',dpi=300);fig.savefig(F/'10_program_cross_cohort.pdf');plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(12,3.7))
 for ax,label in zip(axs,['GSE166466','GSE147009','GSE51993']):
  c,e,s=cohorts[label];levels=['Normal','HSIL','Cancer'] if label=='GSE166466' else ['HPV16_normal','HSIL','Cancer'] if label=='GSE147009' else ['Normal','CIN1','HSIL']
  for i,lev in enumerate(levels):
   vals=s.loc[c.index[c.stage.eq(lev)],'Joint_program'];xx=np.linspace(i-.12,i+.12,len(vals));ax.scatter(xx,vals,s=26,color=['#637f91','#c7a55b','#ad6a50'][i])
  ax.set(xticks=range(3),xticklabels=levels,title=label,ylabel='Fixed program signed mean fractional rank');ax.tick_params(axis='x',rotation=20)
 fig.tight_layout();fig.savefig(F/'11_program_stage_context.png',dpi=300);fig.savefig(F/'11_program_stage_context.pdf');plt.close(fig)
 print('Program stage contrasts:\n',staged[staged.target.eq('Joint_program')].to_string(index=False));print('\nComposition-adjusted program:\n',composition[composition.gene.eq('Joint_program')].to_string(index=False));print('\nComposition diagnostics:\n',pd.DataFrame(diag).to_string(index=False));print('\nProspective program:\n',clinical[clinical.gene.eq('Joint_program')].to_string(index=False))

if __name__=='__main__':main()

import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,P,R,design,ols_all,dump
from extension_utils import prepare_new_bulk,frozen_targets,fit_targets,program_fit,extra_scores,read_meta

def main():
 cohorts=prepare_new_bulk();old,meta,directions=frozen_targets();union=list(directions)
 rows=[];prog=[]
 for label in ['GSE166466','GSE149763','GSE75132']:
  c,e,s=cohorts[label]
  for adj in [True,False]:
   rows+=fit_targets(label,c,e,s,union,directions,adj);prog.append(program_fit(label,c,s,adj))
 result=pd.DataFrame(rows);program=pd.DataFrame(prog)
 result['holm_union_within_cohort_model']=result.groupby(['cohort','model']).p.transform(lambda x:multipletests(x,method='holm')[1])
 result['holm_old_three_within_cohort_model']=np.nan;result['holm_meta_ten_within_cohort_model']=np.nan
 for (_,part) in result.groupby(['cohort','model']):
  for genes,column in [(old,'holm_old_three_within_cohort_model'),(meta,'holm_meta_ten_within_cohort_model')]:
   ix=part.index[part.gene.isin(genes)];result.loc[ix,column]=multipletests(result.loc[ix,'p'],method='holm')[1]
 primary=result.model.eq('proliferation_adjusted')&result.cohort.isin(['GSE166466','GSE149763']);ii=result.index[primary]
 result['holm_union_across_two_new_studies']=np.nan;result.loc[ii,'holm_union_across_two_new_studies']=multipletests(result.loc[ii,'p'],method='holm')[1]
 program['holm_program_across_two_new_studies']=np.nan;jj=program.index[program.model.eq('proliferation_adjusted')&program.cohort.isin(['GSE166466','GSE149763'])]
 program.loc[jj,'holm_program_across_two_new_studies']=multipletests(program.loc[jj,'p'],method='holm')[1]
 allp=pd.concat([result.loc[ii,['cohort','gene','p']],program.loc[jj,['cohort','gene','p']]])
 allp['holm_all_genes_and_program_two_new_studies']=multipletests(allp.p,method='holm')[1]
 result=result.merge(allp[['cohort','gene','holm_all_genes_and_program_two_new_studies']],on=['cohort','gene'],how='left')
 program=program.merge(allp[['cohort','gene','holm_all_genes_and_program_two_new_studies']],on=['cohort','gene'],how='left')
 # Merged conservative P values refer to primary adjusted models only.
 for dd in [result,program]:dd.loc[dd.model.ne('proliferation_adjusted'),'holm_all_genes_and_program_two_new_studies']=np.nan
 result.to_csv(R/'additional_cohort_all_frozen_gene_tests.tsv',sep='\t',index=False);program.to_csv(R/'joint_program_new_cohort_tests.tsv',sep='\t',index=False)
 # Author RPKM sensitivity does not replace the primary read-derived CPM fit.
 c,e,s=cohorts['GSE149763'];er=pd.read_csv(P/'GSE149763_author_log_RPKM.tsv.gz',sep='\t',index_col=0);sr,_=extra_scores(er)
 sens=pd.DataFrame(fit_targets('GSE149763_author_RPKM',c,er,sr,union,directions,True)+[program_fit('GSE149763_author_RPKM',c,sr,True)]);sens['holm_union_plus_program']=multipletests(sens.p,method='holm')[1];sens.to_csv(R/'GSE149763_author_RPKM_sensitivity.tsv',sep='\t',index=False)
 # Platform audit only, no invalid HSIL/SCC regression.
 audit=[]
 for platform in ['GPL570','GPL571']:
  c=read_meta(f'GSE27678-{platform}_series_matrix.txt.gz');c['stage']=np.select([c.title.str.startswith('HSIL'),c.title.str.startswith('LSIL'),c.title.str.startswith('SCC'),c.title.str.startswith('CellLine')],['HSIL','LSIL','SCC','CellLine'],default='Normal');c['platform']=platform;audit.append(c)
 audit=pd.concat(audit);audit.to_csv(P/'GSE27678_platform_pathology_audit.tsv',sep='\t');pd.crosstab(audit.stage,audit.platform).to_csv(R/'GSE27678_platform_pathology_crosstab.tsv',sep='\t')
 for label,(c,e,s) in cohorts.items():print(label,c.stage.value_counts().to_dict(),e.shape,flush=True)
 print('\nOriginal three new-study adjusted results:\n',result[result.gene.isin(old)&result.model.eq('proliferation_adjusted')].to_string(index=False));print('\nJoint program:\n',program.to_string(index=False));print('\nNew meta targets:\n',result[result.gene.isin(meta)&result.model.eq('proliferation_adjusted')][['cohort','gene','beta','p','holm_meta_ten_within_cohort_model','holm_union_across_two_new_studies','direction_concordant']].to_string(index=False))

if __name__=='__main__':main()

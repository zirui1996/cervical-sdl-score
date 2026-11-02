"""Same patients, independent processing resource: explicitly sensitivity only."""
import io,json,zipfile
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import P,R,RAW,dump
from extension_utils import extra_scores,frozen_targets,fit_targets,program_fit

old,meta,directions=frozen_targets();rows=[];audits=[]
for label in ['GSE223804','GSE149763']:
 if label=='GSE223804':
  c=pd.read_csv(P/'GSE223804_metadata.tsv',sep='\t',index_col=0);c=c.set_index('sample',drop=False)
 else:c=pd.read_csv(P/f'{label}_extension_metadata.tsv',sep='\t',index_col=0)
 with zipfile.ZipFile(RAW/(label+'_ARCHS4_expression.zip')) as z:
  assert z.testzip() is None
  counts=pd.read_csv(io.BytesIO(z.read('matrix.tsv')),sep='\t',index_col=0)
 assert set(counts.columns)==set(c.index) and (counts.values>=0).all() and np.isfinite(counts.values).all()
 audits.append({'cohort':label,'samples':len(c),'all_sample_ids_match':True,'duplicate_symbol_rows':int(counts.index.duplicated().sum()),'is_new_patient_cohort':False})
 for agg in ['sum','max']:
  ct=counts.groupby(level=0).agg(agg);e=np.log2(ct.div(ct.sum(axis=0),axis=1)*1e6+1);s,cov=extra_scores(e)
  ep=P/f'{label}_ARCHS4_{agg}_logCPM.tsv.gz';e.to_csv(ep,sep='\t',compression='gzip');s.to_csv(P/f'{label}_ARCHS4_{agg}_scores.tsv',sep='\t')
  recs=fit_targets(label,c,e,s,list(directions),directions,True)+[program_fit(label,c,s,True)]
  for rr in recs:rr['processing']='ARCHS4_'+agg;rr['role']='same_patient_processing_sensitivity'
  rows+=recs
result=pd.DataFrame(rows);result['holm_targets_and_program_within_cohort_processing']=result.groupby(['cohort','processing']).p.transform(lambda x:multipletests(x,method='holm')[1]);result.to_csv(R/'uniform_RNA_reprocessing_sensitivity.tsv',sep='\t',index=False)
primary=result[result.processing.eq('ARCHS4_sum')].copy();third=pd.read_csv(R/'new_three_study_full_multiplicity.tsv',sep='\t');third=third[third.cohort.eq('GSE147009')].copy();third['processing']='ARCHS4_sum';third['role']='new_held_out_tissue'
uniform=pd.concat([primary,third],ignore_index=True);uniform['holm_all_targets_and_program_uniform_RNA']=multipletests(uniform.p,method='holm')[1];uniform.to_csv(R/'uniform_RNA_three_study_all_targets.tsv',sep='\t',index=False)
dump('uniform_RNA_source_audit.json',{'cohort_audits':audits,'reference_cohort':'GSE147009','mapping':'Repeated ARCHS4 gene symbols: sum primary sensitivity convention; max separate sensitivity','role':'GSE223804 and GSE149763 are reprocessings, not additional patients/confirmatory studies','processing_service':json.loads((RAW/'ARCHS4_status.json').read_text())})
print(uniform[uniform.gene.eq('Joint_program')][['cohort','beta','ci_low','ci_high','p','holm_all_targets_and_program_uniform_RNA','direction_concordant']].to_string(index=False));print('\nAll original three genes, uniform RNA resource:\n',uniform[uniform.gene.isin(old)][['cohort','gene','beta','p','direction_concordant']].to_string(index=False))

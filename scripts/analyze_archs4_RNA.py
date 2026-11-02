from pathlib import Path
import json,zipfile,io
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import P,R,RAW,dump
from extension_utils import extra_scores,frozen_targets,fit_targets,program_fit

with zipfile.ZipFile(RAW/'GSE147009_ARCHS4_expression.zip') as z:
 assert z.testzip() is None
 counts=pd.read_csv(io.BytesIO(z.read('matrix.tsv')),sep='\t',index_col=0)
expected=json.loads((RAW/'GSE147009_ARCHS4_metadata_search.json').read_text())['samples']
assert set(counts.columns)==set(expected) and len(expected)==24
assert (counts.to_numpy()>=0).all() and np.isfinite(counts.to_numpy()).all()
source_rows=len(counts);duplicate_rows=int(counts.index.duplicated().sum())
duplicate_audit=counts.loc[counts.index.duplicated(keep=False)].copy()
audit_path=P/'GSE147009_ARCHS4_duplicate_symbol_rows.tsv.gz';audit_tmp=P/'GSE147009_ARCHS4_duplicate_symbol_rows.atomic.tsv.gz'
duplicate_audit.to_csv(audit_tmp,sep='\t',compression='gzip');assert pd.read_csv(audit_tmp,sep='\t',index_col=0).shape==duplicate_audit.shape;audit_tmp.replace(audit_path)
# The public extraction includes repeated reference gene symbols. Retain the
# source, sum their distinct feature counts, and explicitly test max aggregation
# as a sensitivity because source Ensembl identifiers were not returned.
max_counts=counts.groupby(level=0).max()
counts=counts.groupby(level=0).sum()
meta=pd.DataFrame(index=expected);meta.index.name='sample';access=meta.index.str.extract(r'GSM(\d+)',expand=False).astype(int)
meta['stage']=np.select([access<=4412232,access<=4412238,access<=4412244],['HPV_negative_normal','HPV16_normal','HSIL'],default='Cancer')
# Exact ranges are copied from the official GEO series sample list.
meta['group']=meta.stage.eq('Cancer').astype(int);meta['primary_included']=meta.stage.isin(['HSIL','Cancer'])
assert meta.stage.value_counts().eq(6).all()
libraries=counts.sum(axis=0);e=np.log2(counts.div(libraries,axis=1)*1e6+1);s,cov=extra_scores(e)
meta.to_csv(P/'GSE147009_extension_metadata.tsv',sep='\t');e.to_csv(P/'GSE147009_extension_gene_log_expression.tsv.gz',sep='\t',compression='gzip');s.to_csv(P/'GSE147009_extension_scores.tsv',sep='\t');cov.to_csv(R/'GSE147009_extension_score_coverage.tsv',sep='\t',index=False)
dump('GSE147009_ARCHS4_scale_audit.json',{'samples':len(meta),'expected_sample_set_complete':True,'source_rows':source_rows,'source_duplicate_symbol_rows':duplicate_rows,'unique_symbols_after_sum':len(counts),'mapping':'sum repeated source symbols; max aggregation sensitivity; source Ensembl identifiers were not returned','library_sums':libraries.to_dict(),'integer_values':bool(np.allclose(counts,counts.round())),'normalization':'log2(1+CPM) from ARCHS4 uniformly processed estimated counts','processing_resource':'ARCHS4; Kallisto-based; not author circRNA matrix','service_status':json.loads((RAW/'ARCHS4_status.json').read_text())})
old,metagenes,directions=frozen_targets();rows=[];program=[]
for adj in [True,False]:
 rows+=fit_targets('GSE147009',meta,e,s,list(directions),directions,adj);program.append(program_fit('GSE147009',meta,s,adj))
d=pd.DataFrame(rows);g=pd.DataFrame(program)
em=np.log2(max_counts.div(max_counts.sum(axis=0),axis=1)*1e6+1);sm,_=extra_scores(em)
maxsens=pd.DataFrame(fit_targets('GSE147009_max_symbol_sensitivity',meta,em,sm,list(directions),directions,True)+[program_fit('GSE147009_max_symbol_sensitivity',meta,sm,True)])
maxsens['holm_union_and_program']=multipletests(maxsens.p,method='holm')[1];maxsens.to_csv(R/'GSE147009_max_duplicate_symbol_sensitivity.tsv',sep='\t',index=False)
d['holm_union_within_cohort_model']=d.groupby('model').p.transform(lambda x:multipletests(x,method='holm')[1])
d['holm_old_three_within_cohort_model']=np.nan;d['holm_meta_ten_within_cohort_model']=np.nan
for model,part in d.groupby('model'):
 for targets,column in [(old,'holm_old_three_within_cohort_model'),(metagenes,'holm_meta_ten_within_cohort_model')]:
  ix=part.index[part.gene.isin(targets)];d.loc[ix,column]=multipletests(d.loc[ix,'p'],method='holm')[1]
d.to_csv(R/'GSE147009_ARCHS4_frozen_gene_tests.tsv',sep='\t',index=False);g.to_csv(R/'GSE147009_ARCHS4_program_tests.tsv',sep='\t',index=False)
prior=pd.read_csv(R/'additional_cohort_all_frozen_gene_tests.tsv',sep='\t');previousprogram=pd.read_csv(R/'joint_program_new_cohort_tests.tsv',sep='\t')
allgenes=pd.concat([prior[prior.cohort.isin(['GSE166466','GSE149763'])],d],ignore_index=True)
allprogram=pd.concat([previousprogram[previousprogram.cohort.isin(['GSE166466','GSE149763'])],g],ignore_index=True)
allgenes=allgenes[allgenes.model.eq('proliferation_adjusted')].copy();allprogram=allprogram[allprogram.model.eq('proliferation_adjusted')].copy()
allprogram['holm_program_across_three_new_studies']=multipletests(allprogram.p,method='holm')[1]
allgenes['holm_union_genes_across_three_new_studies']=multipletests(allgenes.p,method='holm')[1]
combined=pd.concat([allgenes,allprogram],ignore_index=True);combined['holm_all_genes_and_program_three_new_studies']=multipletests(combined.p,method='holm')[1]
combined.to_csv(R/'new_three_study_full_multiplicity.tsv',sep='\t',index=False)
print(combined[['cohort','gene','beta','ci_low','ci_high','p','holm_all_genes_and_program_three_new_studies','direction_concordant']].to_string(index=False));print('\nPrograms across three new studies:\n',allprogram.to_string(index=False))

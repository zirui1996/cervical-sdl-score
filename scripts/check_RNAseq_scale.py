"""Document ambiguous author-declared CPM scale; retain primary and labeled sensitivity."""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import RAW,P,R,design,ols_all,dump
import json
raw=pd.read_csv(RAW/'GSE223804_cpm_CIN_CR_Norm.txt.gz',sep='\t',index_col=0)
c=pd.read_csv(P/'GSE223804_metadata.tsv',sep='\t',index_col=0);c=c.loc[c.primary_included]
s=pd.read_csv(P/'GSE223804_proliferation_scores.tsv',sep='\t',index_col=0);genes=json.loads((R/'RNAseq_candidate_lock.json').read_text())['genes']
# Column scaling preserves within-sample ranks, so the frozen proliferation scores are unchanged.
normalized=np.log2(raw.div(raw.sum(axis=0),axis=1)*1e6+1)
res=ols_all(normalized.loc[genes],design(c,s));res['holm_p']=multipletests(res.p,method='holm')[1]
res['analysis']='sensitivity_only_column_total_rescaling_to_one_million_then_log2'
res.to_csv(R/'RNAseq_column_scale_sensitivity.tsv',sep='\t',index=False)
dump('RNAseq_scale_audit.json',{'author_declared_units':'CPM','all_values_integer':bool(np.equal(raw,np.floor(raw)).all().all()),'column_sum_min':float(raw.sum().min()),'column_sum_max':float(raw.sum().max()),'plain_full_matrix_CPM_sum_expectation':1000000,'interpretation':'File provenance declares CPM, but integer values and widely varying column totals make the precise normalization uncertain. Primary results preserve the supplied values. Column-total scaling is only a sensitivity and cannot establish that the input was raw counts or a complete transcriptome. The RNA cohort is not decisive validation without clarification.'})
print(res.to_string(index=False))

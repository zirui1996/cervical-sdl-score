"""Scientific checks for frozen choices, source roles and reported inference."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,P,R,design,gene_sets
from extension_utils import frozen_targets,extra_scores

checks=[]
def check(name,value):
 assert bool(value),name
 checks.append({'check':name,'passed':True})

lock=json.loads((R/'meta_discovery_candidate_lock.json').read_text());prog=json.loads((R/'joint_program_lock.json').read_text());d=pd.read_csv(R/'meta_discovery_all_genes.tsv',sep='\t')
check('Joint discovery tests all 6388 fixed common measured genes',len(d)==6388 and not d.gene.duplicated().any())
check('Joint discovery BH covers the whole tested universe',np.allclose(d.combined_fdr,multipletests(d.combined_p,method='fdr_bh')[1]))
source=pd.read_csv(R/'genomewide_unadjusted_adjusted_coefficients.tsv.gz',sep='\t');labels=['GSE63514','GSE138080','GSE7803'];ns=np.array([60,25,28])
old=source.pivot(index='gene',columns='cohort',values='p_adjusted').reindex(d.gene)[labels];beta=source.pivot(index='gene',columns='cohort',values='beta_adjusted').reindex(d.gene)[labels]
zz=np.sum(np.sign(beta.to_numpy())*norm.isf(np.clip(old.to_numpy(),1e-300,1)/2)*np.sqrt(ns),axis=1)/np.sqrt(ns.sum())
check('Signed Stouffer statistic independently reconstructed',np.allclose(d.combined_z,zz))
chosen=d[d.same_direction_all_three&d.combined_fdr.lt(.05)&d.mean_absolute_standardized_beta.ge(.5)].sort_values(['combined_p','gene']).head(10)
check('Frozen ten genes follow the prewritten joint rule',[x['gene'] for x in lock['candidates']]==chosen.gene.tolist())
check('Meta ranking plan hash preserved',hashlib.sha256((ROOT/'META_DISCOVERY_EXTENSION_v1.md').read_bytes()).hexdigest()==lock['plan_sha256'])
check('Program definition plan hash preserved',hashlib.sha256((ROOT/'JOINT_PROGRAM_TEST_v1.md').read_bytes()).hexdigest()==prog['plan_sha256'])
check('Program uses exactly ten frozen equal-direction targets',{x['gene'] for x in prog['candidates']}==set(chosen.gene) and all(x['direction']==-1 for x in prog['candidates']))
check('Targets do not overlap frozen proliferation-score genes',not(set(chosen.gene)&set.union(*gene_sets().values())))
alltests=pd.read_csv(R/'new_three_study_full_multiplicity.tsv',sep='\t')
check('New combined testing includes all 42 hypotheses',len(alltests)==42 and not alltests[['cohort','gene']].duplicated().any())
check('All-new-family Holm independently recalculated',np.allclose(multipletests(alltests.p,method='holm')[1],alltests.holm_all_genes_and_program_three_new_studies))
for label,n,counts in [('GSE166466',13,{'Cancer':7,'HSIL':6}),('GSE149763',6,{'HSIL':3,'Cancer':3}),('GSE147009',12,{'HSIL':6,'Cancer':6})]:
 c=pd.read_csv(P/f'{label}_extension_metadata.tsv',sep='\t',index_col=0);cc=c.loc[c.primary_included];e=pd.read_csv(P/f'{label}_extension_gene_log_expression.tsv.gz',sep='\t',index_col=0);s,cov=extra_scores(e)
 check(label+' primary pathology and counts',len(cc)==n and cc.stage.value_counts().to_dict()==counts)
 check(label+' measures every frozen program symbol',int(cov.loc[cov.score.eq('Joint_program'),'available_genes'].iloc[0])==10)
 v=s.loc[cc.index,'Joint_program'];y=(v-v.mean())/v.std(ddof=1);fit=sm.OLS(y,design(cc,s)).fit(cov_type='HC3',use_t=True);rr=alltests[alltests.cohort.eq(label)&alltests.gene.eq('Joint_program')].iloc[0]
 check(label+' program beta and P match independent statsmodels',np.isclose(fit.params['group'],rr.beta) and np.isclose(fit.pvalues['group'],rr.p))
 for gene in ['FCGBP','GOT2']:
  fit=sm.OLS(e.loc[gene,cc.index],design(cc,s)).fit(cov_type='HC3',use_t=True);rr=alltests[alltests.cohort.eq(label)&alltests.gene.eq(gene)].iloc[0]
  check(label+' '+gene+' HC3 independently reproduced',np.isclose(fit.params['group'],rr.beta) and np.isclose(fit.pvalues['group'],rr.p))
stage=pd.read_csv(R/'stage_context_all_contrasts.tsv',sep='\t')
check('Stage comparison multiplicity includes all 252 planned target contrasts',len(stage)==252 and np.allclose(multipletests(stage.p,method='fdr_bh')[1],stage.extension_stage_fdr))
q=pd.read_csv(R/'GSE208653_QC_summary.tsv',sep='\t');pb=pd.read_csv(R/'GSE208653_specimen_compartment_programs.tsv',sep='\t');states=pd.read_csv(R/'GSE208653_epithelial_marker_states.tsv',sep='\t')
check('Single-cell QC accounts for all 83959 source cells',len(q)==9 and q.input_cells.sum()==83959 and np.all(q.input_cells==q.retained_cells+q.excluded_cells))
check('76733 retained cells are consistently summarized',q.retained_cells.sum()==76733 and pb.loc[pb.compartment.eq('All cells'),'cells'].sum()==76733)
check('ADC remains a separate ninth specimen',q.stage.value_counts().to_dict()=={'HPV_negative_normal':2,'HPV_positive_normal':2,'HSIL':2,'SCC':2,'ADC':1})
check('Marker compartments partition QC cells per donor',pb[~pb.compartment.eq('All cells')].groupby('donor').cells.sum().equals(pb[pb.compartment.eq('All cells')].set_index('donor').cells.sort_index()))
ep=pb[pb.compartment.eq('Epithelial')].set_index('donor').cells.sort_index();check('Epithelial marker states partition epithelial cells',states.groupby('donor').cells.sum().sort_index().equals(ep))
check('Sub-30-cell state scores are not promoted',states.loc[states.cells.lt(30),'Joint_program'].isna().all())
clinical=pd.read_csv(R/'prospective_endpoint_transfer.tsv',sep='\t');c=pd.read_csv(P/'GSE75132_extension_metadata.tsv',sep='\t');check('Prospective counts use mutually consistent source codes',c.category_code.value_counts().to_dict()=={2:20,0:11,1:10})
check('Prospective candidate/program family corrected together',len(clinical)==14 and np.allclose(multipletests(clinical.p,method='holm')[1],clinical.holm_union_and_program))
audit=json.loads((R/'GSE75132_label_audit.json').read_text());endpoint=pd.read_csv(R/'prospective_endpoint_definition_sensitivity.tsv',sep='\t')
check('Mixed-grade endpoint caveat retained',audit['strict_CIN3plus_endpoint_not_established'] and audit['disease_state_in_progression_group']['moderate dysplasia']==4)
check('Restricted annotation subset has 26 samples and does not reassign moderate records',set(endpoint.loc[endpoint.endpoint_definition.str.startswith('Restricted'),'n'])=={26})
check('All 28 endpoint-definition sensitivity tests corrected together',len(endpoint)==28 and np.allclose(multipletests(endpoint.p,method='holm')[1],endpoint.holm_all_28_endpoint_tests))
e=pd.read_csv(P/'GSE75132_gene_log_expression.tsv.gz',sep='\t',index_col=0);ex=pd.read_csv(P/'GSE75132_extension_gene_log_expression.tsv.gz',sep='\t',index_col=0)
check('Ratios below one retained as negative log ratios',e.min().min()<0 and np.allclose(e,ex))
check('Original RNA cache repair preserved three-gene results',json.loads((R/'processed_RNA_cache_repair.json').read_text())['original_three_gene_beta_and_p_reproduced'])
check('Original five candidates still reported as failed',not pd.read_csv(R/'external_validation.tsv',sep='\t').replicated.any())
check('Frozen new RNA three-gene results still unconfirmed',not pd.read_csv(R/'new_RNAseq_validation.tsv',sep='\t').replicated.any())
for f in sorted(P.glob('*.tsv.gz')):
 if '.atomic.' in f.name:continue
 with gzip.open(f,'rb') as h:
  while h.read(1024*1024):pass
 check('Processed gzip reopens completely '+f.name,True)
(R/'extension_verification.json').write_text(json.dumps({'all_passed':True,'checks':checks},indent=2))
print('PASS',len(checks),'extension scientific, multiplicity, labeling and matrix checks')

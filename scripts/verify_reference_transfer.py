"""Independent rule, donor, source and mixture checks; no cell-wise inference."""
import hashlib,json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,P,R,design
from extension_utils import common_universe,extra_scores
from analyze_stage_and_composition import load_all

checks=[]
def check(name,value):
    assert bool(value),name
    checks.append({'check':name,'passed':True})

nested=pd.read_csv(R/'nested_discovery_rule_summary.tsv',sep='\t')
membership=pd.read_csv(R/'nested_discovery_rule_membership.tsv',sep='\t')
ranking=pd.read_csv(R/'nested_discovery_rule_all_rankings.tsv.gz',sep='\t')
cohorts=load_all()
for _,row in nested.iterrows():
    d=ranking[ranking.heldout.eq(row.heldout)]
    check(row.heldout+' full-universe BH',len(d)==6388 and np.allclose(multipletests(d.p,method='fdr_bh')[1],d.fdr))
    selected=d[d.same_sign&d.fdr.lt(.05)&d.mean_abs_beta_standardized.ge(.5)].sort_values(['p','gene']).head(10)
    check(row.heldout+' selection restricted to two training-study statistics',selected.gene.tolist()==row.selected_genes.split(';') and row.heldout not in row.training.split(';'))
    c,e,s=cohorts[row.heldout];cc=c.loc[c.primary_included];ranks=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True)
    v=ranks.loc[selected.gene].mul(selected.set_index('gene').direction,axis=0).mean(axis=0);v=(v-v.mean())/v.std(ddof=1)
    fit=sm.OLS(v,design(cc,s)).fit(cov_type='HC3',use_t=True)
    check(row.heldout+' held-out score independently recomputed',np.isclose(fit.params['group'],row.beta) and np.isclose(fit.pvalues['group'],row.p,rtol=1e-8))
check('All three stress-test folds retained and jointly corrected',len(nested)==3 and np.allclose(multipletests(nested.p,method='holm')[1],nested.holm_three_fold_tests))
check('Original method-transfer plan unchanged',json.loads((R/'nested_discovery_rule_record.json').read_text())['plan_sha256']==hashlib.sha256((ROOT/'REFERENCE_METHOD_TRANSFER_v1.md').read_bytes()).hexdigest())
mix=pd.read_csv(R/'captured_RNA_mixture_pair_decomposition.tsv',sep='\t');coverage=pd.read_csv(R/'captured_RNA_mixture_coverage.tsv',sep='\t')
check('All four descriptive crossing pairs retained without independence claim',len(mix)==4 and mix.HSIL_specimen.nunique()==2 and mix.SCC_specimen.nunique()==2)
check('Nonlinear-score decomposition adds exactly to observed contrast',np.allclose(mix.captured_RNA_composition_component+mix.within_compartment_profile_component,mix.observed_restricted_difference,atol=1e-14))
check('Restricted-mixture reconstruction and coverage recorded',coverage.reconstruction_max_absolute_CPM_error.max()<1e-8 and coverage.retained_cell_fraction.min()>.95 and coverage.retained_captured_UMI_fraction.min()>.95)
qc=pd.read_csv(R/'EMTAB12305_QC_summary.tsv',sep='\t');pseudobulk=pd.read_csv(P/'EMTAB12305_compartment_pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
summary=pd.read_csv(R/'EMTAB12305_specimen_compartment_programs.tsv',sep='\t')
check('Nine cervical specimens; metastatic node excluded',len(qc)==9 and qc.stage.value_counts().to_dict()=={'SCC':4,'Adjacent_normal':3,'HSIL':2} and not qc.specimen.eq('L1').any())
check('New QC accounts for all input and retained cells',qc.input_cells.sum()==106998 and qc.retained_cells.sum()==93556 and (qc.retained_cells+qc.excluded_cells==qc.input_cells).all())
check('Native matrices checked for nonnegative integer units',qc.nonnegative_integer_counts_verified.all())
for specimen in qc.specimen:
    sub=summary[summary.specimen.eq(specimen)]
    check(specimen+' compartment counts partition all retained cells',sub.loc[~sub.compartment.eq('All cells'),'cells'].sum()==sub.loc[sub.compartment.eq('All cells'),'cells'].iloc[0])
    check(specimen+' pseudobulk UMI counts partition exactly',np.array_equal(pseudobulk[[x for x in pseudobulk if x.startswith(specimen+'|') and x!=specimen+'|All cells']].sum(axis=1).to_numpy(),pseudobulk[specimen+'|All cells'].to_numpy()))
e=np.log2(pseudobulk.div(pseudobulk.sum(),axis=1)*1e6+1);scores,cov=extra_scores(e)
good=summary.cells.ge(30)
check('New frozen scores independently reconstruct from count sums',np.allclose(summary.loc[good,'Joint_program'],scores.loc[summary.loc[good,'key'],'Joint_program']))
check('Below-thirty-cell compartments cannot support a score',summary.loc[~good,'Joint_program'].isna().all())
epi=summary[summary.compartment.eq('Epithelial')&summary.stage.isin(['HSIL','SCC'])]
check('All six primary specimens have adequate epithelial coverage',len(epi)==6 and epi.cells.ge(30).all())
check('New source directional result is preserved without a cell-wise P',epi.loc[epi.stage.eq('SCC'),'Joint_program'].min()>epi.loc[epi.stage.eq('HSIL'),'Joint_program'].max())
record=json.loads((R/'EMTAB12305_analysis_summary.json').read_text())
check('Mechanism, malignant identity and patient-pairing boundaries retained',record['no_cellwise_inference'] and record['not_malignant_cell_identity'] and record['patient_pairing_not_inferred_from_equal_ages'])
if (R/'EMTAB12305_epithelial_marker_states.tsv').exists():
    states=pd.read_csv(R/'EMTAB12305_epithelial_marker_states.tsv',sep='\t')
    check('Marker states partition epithelial cells in all primary specimens',states.groupby('specimen').cells.sum().sort_index().equals(epi.set_index('specimen').cells.sort_index()))
    check('Small epithelial states retain missing scores',states.loc[states.cells.lt(30),'Joint_program'].isna().all())
    check('State follow-up is explicitly post hoc',json.loads((R/'EMTAB12305_marker_states_record.json').read_text())['posthoc_after_new_program_results'])
(R/'reference_transfer_verification.json').write_text(json.dumps({'passed':True,'checks':checks},indent=2))
print(len(checks),'method-transfer/new-source scientific checks passed')

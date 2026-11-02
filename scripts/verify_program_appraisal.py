"""Independent scientific checks for reference panels and uncertainty estimates."""
import json,hashlib
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,R
from extension_utils import frozen_targets,common_universe
from analyze_stage_and_composition import load_all
from appraise_program_distinctiveness import primary,stage

checks=[]
def check(name,condition):
    assert bool(condition),name
    checks.append({'check':name,'passed':True})

members=pd.read_csv(R/'program_reference_panel_memberships.tsv.gz',sep='\t')
fits=pd.read_csv(R/'program_reference_all_fits.tsv.gz',sep='\t')
summary=pd.read_csv(R/'program_distinctiveness_reference_summary.tsv',sep='\t')
old,targets,directions=frozen_targets();excluded=set(old+targets)
check('Exactly 6000 panels before four context evaluations',len(members)==6000 and members.groupby('family').size().eq(2000).all())
check('Every reference panel has ten unique non-target genes',members.genes.map(lambda s:len(set(s.split(';')))==10 and not(set(s.split(';'))&excluded)).all())
check('All 24000 attempted evaluations retained including unscorable draws',len(fits)==24000 and len(summary)==12)
check('Sub-eight-gene panels cannot be scorable',not fits.loc[fits.measured_genes.lt(8),'valid'].any())
for _,row in summary.iterrows():
    d=fits[(fits.cohort==row.cohort)&(fits.context==row.context)&(fits.family==row.family)]
    good=d.loc[d.valid];tail=(1+(good.t>=row.fixed_t).sum())/(1+len(good))
    check(row.cohort+' '+row.context+' '+row.family+' reference summary',len(d)==2000 and len(good)==row.scorable_draws and np.isclose(tail,row.add_one_reference_tail))
check('Twelve descriptive reference tails corrected without dropping families',np.allclose(multipletests(summary.add_one_reference_tail,method='holm')[1],summary.holm_12_descriptive_reference_tails))
cohorts=load_all()
for label,context in [('GSE166466','invasion'),('GSE149763','invasion'),('GSE147009','invasion'),('GSE166466','late_minus_early')]:
    c,e,s=cohorts[label];cc,X,contrast=primary(c,s) if context=='invasion' else stage(c,s)
    rank=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True)
    for family in members.family.unique():
        row=members[(members.family==family)&members.draw.eq(0)].iloc[0]
        genes=[g for g in row.genes.split(';') if g in rank.index]
        y=-rank.loc[genes].mean(axis=0);y=(y-y.mean())/y.std(ddof=1)
        fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True).t_test(contrast)
        ref=fits[(fits.cohort==label)&(fits.context==context)&(fits.family==family)&fits.draw.eq(0)].iloc[0]
        check(label+' '+context+' '+family+' independent reference score/HC3',np.isclose(float(np.asarray(fit.effect).item()),ref.beta) and np.isclose(float(np.asarray(fit.pvalue).item()),ref.p,rtol=1e-8,atol=1e-12))
rna=pd.read_csv(R/'RNA_program_compatible_effect_ranges.tsv',sep='\t')
check('RNA intervals cross zero and include Clariom point; no absence claim',len(rna)==3 and rna.interval_contains_zero.all() and rna.interval_contains_Clariom_point.all())
pooled=pd.read_csv(R/'RNA_program_exploratory_pooled_estimation.tsv',sep='\t')
for _,row in pooled.iterrows():
    d=rna[rna.cohort.isin(row.cohorts.split(';'))];weights=1/(d.se_HC3**2+row.tau_squared_PM)
    beta=np.average(d.beta,weights=weights)
    check(row.context+' inverse-variance mean and conservative KH interval',np.isclose(beta,row.beta) and row.modified_KH_se>=np.sqrt(1/weights.sum())-1e-12 and np.isclose(row.ci_high-row.beta,stats.t.ppf(.975,len(d)-1)*row.modified_KH_se))
record=json.loads((R/'program_distinctiveness_run_record.json').read_text())
check('Written plan and original meta inputs are preserved',record['plan_sha256']==hashlib.sha256((ROOT/'PROGRAM_DISTINCTIVENESS_v1.md').read_bytes()).hexdigest() and record['meta_input_sha256']==hashlib.sha256((R/'meta_discovery_all_genes.tsv').read_bytes()).hexdigest())
check('Post hoc timing is explicit',record['posthoc_after_original_outcomes_known'])
matching=pd.read_csv(R/'program_reference_realized_matching.tsv.gz',sep='\t')
features=[k for k in matching if k.startswith('scaled_difference_')]
quality=matching[features].agg(['mean','std','min','max']).T.reset_index(names='feature')
quality['acceptable_absolute_mean_below_0p25']=quality['mean'].abs()<.25
quality.to_csv(R/'program_reference_matching_quality.tsv',sep='\t',index=False)
check('Evidence-strength mismatch flagged rather than called exchangeable',quality.loc[quality.feature.eq('scaled_difference_abs_combined_z'),'mean'].iloc[0]<-1)
(R/'program_appraisal_verification.json').write_text(json.dumps({'passed':True,'checks':checks,'matching_caveat':'The top ten selected genes have extreme discovery evidence; nearest-neighbour controls remain weaker, mean combined-Z difference below -1 SD. These are not evidence-equivalent controls.'},indent=2))
print(len(checks),'appraisal scientific checks passed; matching limitations retained')

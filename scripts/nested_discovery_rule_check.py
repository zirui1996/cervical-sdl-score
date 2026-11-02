"""Reselect in two old studies; stress-test rule in the third, post hoc."""
import json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,R,design,ols_all
from extension_utils import common_universe,frozen_targets
from analyze_stage_and_composition import load_all

def main():
    cohorts=load_all();labels=['GSE63514','GSE138080','GSE7803']
    ns={'GSE63514':60,'GSE138080':25,'GSE7803':28}
    allstats=pd.read_csv(R/'genomewide_unadjusted_adjusted_coefficients.tsv.gz',sep='\t')
    genes=pd.read_csv(R/'genomewide_common_measured_universe.tsv',sep='\t').gene.tolist()
    old,original,direction=frozen_targets();summary=[];membership=[];ranks=[]
    for heldout in labels:
        train=[s for s in labels if s!=heldout]
        dd=allstats[allstats.cohort.isin(train)]
        b=dd.pivot(index='gene',columns='cohort',values='beta_adjusted').loc[genes,train]
        p=dd.pivot(index='gene',columns='cohort',values='p_adjusted').loc[genes,train]
        sb=dd.pivot(index='gene',columns='cohort',values='standardized_beta_adjusted').loc[genes,train]
        w=np.sqrt([ns[x] for x in train]);z=(np.sign(b.to_numpy())*stats.norm.isf(np.clip(p.to_numpy(),1e-300,1)/2))@w/np.linalg.norm(w)
        pp=2*stats.norm.sf(abs(z))
        table=pd.DataFrame({'gene':genes,'combined_z':z,'p':pp,'fdr':multipletests(pp,method='fdr_bh')[1],'same_sign':np.sign(b).nunique(axis=1).eq(1).to_numpy(),'mean_abs_beta_standardized':sb.abs().mean(axis=1).to_numpy(),'direction':np.sign(b.iloc[:,0]).astype(int).to_numpy()})
        table['heldout']=heldout;ranks.append(table)
        qualified=table[table.same_sign&table.fdr.lt(.05)&table.mean_abs_beta_standardized.ge(.5)].sort_values(['p','gene'])
        selected=qualified.head(10);rec={'heldout':heldout,'training':';'.join(train),'qualified_genes':len(qualified),'n_test':ns[heldout],'selected_genes':';'.join(selected.gene),'overlap_with_original_ten':len(set(selected.gene)&set(original)),'status':'not_evaluable_less_than_ten' if len(selected)<10 else 'tested'}
        for i,row in selected.reset_index(drop=True).iterrows():membership.append({'heldout':heldout,'rank':i+1,**row.to_dict()})
        if len(selected)==10:
            c,e,s=cohorts[heldout];cc=c.loc[c.primary_included]
            rank=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True)
            have=[g for g in selected.gene if g in rank.index];rec['measured']=len(have)
            assert len(have)==10
            signs=selected.set_index('gene').direction
            score=rank.loc[have].mul(signs.loc[have],axis=0).mean(axis=0)
            y=pd.DataFrame([(score-score.mean())/score.std(ddof=1)],index=['Fold_specific_score'])
            fit=ols_all(y,design(cc,s)).iloc[0].to_dict();fit.pop('fdr');rec.update(fit)
        else:rec['p']=1.
        summary.append(rec)
    out=pd.DataFrame(summary);out['holm_three_fold_tests']=multipletests(out.p,method='holm')[1]
    out['analysis_role']='posthoc_discovery_rule_stress_test; fold-specific sets; not confirmation of original ten genes'
    out.to_csv(R/'nested_discovery_rule_summary.tsv',sep='\t',index=False)
    pd.DataFrame(membership).to_csv(R/'nested_discovery_rule_membership.tsv',sep='\t',index=False)
    pd.concat(ranks).to_csv(R/'nested_discovery_rule_all_rankings.tsv.gz',sep='\t',index=False,compression='gzip')
    (R/'nested_discovery_rule_record.json').write_text(json.dumps({'training_outcomes_only_in_each_reselection':True,'original_signature_changed':False,'test_data_previously_seen_in_original_project':True,'frozen_6388_universe_uses_original_outcome_independent_expression_filters':True,'not_a_blinded_or_fully_raw_nested_validation':True,'plan_sha256':hashlib.sha256((ROOT/'REFERENCE_METHOD_TRANSFER_v1.md').read_bytes()).hexdigest()},indent=2))
    print(out.to_string(index=False),flush=True)

if __name__=='__main__':main()

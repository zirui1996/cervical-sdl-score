"""Post hoc fixed-program reference comparisons; never selects a replacement."""
import hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats,optimize
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import ROOT,R,F,design
from extension_utils import common_universe,frozen_targets
from analyze_stage_and_composition import load_all

SEED=2026102907
B=2000

def fit_many(Y,X,contrast):
    """Y=sample x panel; normalize outcomes exactly as existing program fits."""
    Y=np.asarray(Y,float)
    sd=np.std(Y,axis=0,ddof=1)
    valid=np.isfinite(Y).all(axis=0)&(sd>1e-12)
    y=(Y[:,valid]-Y[:,valid].mean(axis=0))/sd[valid]
    a=np.asarray(X,float);v=np.asarray(contrast,float)
    inv=np.linalg.inv(a.T@a);coef=inv@a.T@y
    h=np.einsum('ij,jk,ik->i',a,inv,a)
    weights=a@inv@v
    se=np.sqrt(np.sum((weights[:,None]*(y-a@coef)/(1-h[:,None]))**2,axis=0))
    effect=v@coef;tt=effect/se;df=len(a)-a.shape[1]
    out=pd.DataFrame(index=range(Y.shape[1]),columns=['beta','se_HC3','t','p','ci_low','ci_high'],dtype=float)
    crit=stats.t.ppf(.975,df)
    out.loc[valid,:]=np.column_stack([effect,se,tt,2*stats.t.sf(abs(tt),df),effect-crit*se,effect+crit*se])
    out['valid']=valid;out['residual_df']=df
    return out

def primary(c,s):
    cc=c.loc[c.primary_included].copy()
    return cc,design(cc,s),[0,1,0]

def stage(c,s):
    cc=c.loc[c.stage.isin(['Normal','HSIL','Cancer'])].copy()
    v=s.loc[cc.index,'Composite']
    xx=pd.DataFrame({'const':1.,'HSIL':cc.stage.eq('HSIL').astype(float),'Cancer':cc.stage.eq('Cancer').astype(float),'proliferation_z':(v-v.mean())/v.std(ddof=1)},index=cc.index)
    return cc,xx,[0,-2,1,0]

def random_effects(d,label):
    """Paule-Mandel with modified KH factor >= 1 and t_(k-1)."""
    y=d.beta.to_numpy();vv=d.se_HC3.to_numpy()**2;k=len(y)
    def Q(tau):
        w=1/(vv+tau);mu=np.dot(w,y)/w.sum()
        return np.dot(w,(y-mu)**2)
    tau=0.
    if Q(0)>k-1:
        hi=1.
        while Q(hi)>k-1:hi*=2
        tau=optimize.brentq(lambda x:Q(x)-(k-1),0,hi)
    w=1/(vv+tau);mu=np.dot(w,y)/w.sum();q=Q(tau)/(k-1)
    se=np.sqrt(max(1.,q)/w.sum());crit=stats.t.ppf(.975,k-1)
    return {'context':label,'k':k,'cohorts':';'.join(d.cohort),'beta':mu,'tau_squared_PM':tau,'KH_q':q,'modified_KH_se':se,'df':k-1,'ci_low':mu-crit*se,'ci_high':mu+crit*se,'p_exploratory':2*stats.t.sf(abs(mu/se),k-1),'interpretation':'secondary_estimation_not_new_confirmation; cohort-specific SD units; small-k uncertainty'}

def main():
    cohorts=load_all();rng=np.random.default_rng(SEED)
    old,targets,directions=frozen_targets();assert all(directions[g]==-1 for g in targets)
    excluded=set(old+targets)
    meta=pd.read_csv(R/'meta_discovery_all_genes.tsv',sep='\t').set_index('gene').sort_index()
    qualified=meta.same_direction_all_three&meta.combined_fdr.lt(.05)&meta.mean_absolute_standardized_beta.ge(.5)&meta.direction.eq(-1)
    negatives=meta.index[qualified].tolist();pool=[g for g in negatives if g not in excluded]
    allpool=[g for g in meta.index if g not in excluded]
    rank_features=[];sd_features=[]
    for label in ['GSE63514','GSE138080','GSE7803']:
        c,e,s=cohorts[label];cc=c.loc[c.primary_included]
        ranks=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True)
        rank_features.append(ranks.loc[meta.index].mean(axis=1));sd_features.append(ranks.loc[meta.index].std(axis=1,ddof=1))
    features=pd.DataFrame({'abs_combined_z':meta.combined_z.abs(),'mean_abs_standardized_beta':meta.mean_absolute_standardized_beta,'mean_discovery_rank':pd.concat(rank_features,axis=1).mean(axis=1),'mean_discovery_rank_sd':pd.concat(sd_features,axis=1).mean(axis=1)})
    scale=features.loc[pool].std(ddof=1);assert (scale>0).all()
    nearest={};candidate_rows=[]
    for target in targets:
        distance=np.sqrt(((features.loc[pool]-features.loc[target])/scale).pow(2).sum(axis=1))
        near=distance.sort_values(kind='stable').head(20)
        nearest[target]=near.index.tolist()
        for gene,dist in near.items():candidate_rows.append({'target':target,'candidate':gene,'distance':dist,**features.loc[gene].to_dict()})
    features.to_csv(R/'program_reference_discovery_features.tsv',sep='\t')
    pd.DataFrame(candidate_rows).to_csv(R/'program_reference_matching_candidates.tsv',sep='\t',index=False)
    memberships=[];panels={};matching=[]
    for family in ['All_discovery','Qualified_negative','Discovery_matched']:
        panels[family]=[]
        for draw in range(B):
            if family!='Discovery_matched':
                genes=rng.choice(allpool if family=='All_discovery' else pool,10,replace=False).tolist()
            else:
                by_target={};used=set()
                for target in rng.permutation(targets):
                    available=[g for g in nearest[target] if g not in used]
                    assert available
                    gene=str(rng.choice(available));by_target[target]=gene;used.add(gene)
                genes=[by_target[t] for t in targets]
                for j,t in enumerate(targets):
                    delta=(features.loc[genes[j]]-features.loc[t])/scale
                    matching.append({'draw':draw,'target':t,'matched_gene':genes[j],'distance':float(np.linalg.norm(delta)),**{'scaled_difference_'+key:value for key,value in delta.items()}})
            assert len(set(genes))==10 and not (set(genes)&excluded)
            panels[family].append(genes)
            memberships.append({'family':family,'draw':draw,'genes':';'.join(genes)})
    pd.DataFrame(memberships).to_csv(R/'program_reference_panel_memberships.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(matching).to_csv(R/'program_reference_realized_matching.tsv.gz',sep='\t',index=False,compression='gzip')
    panel_fits=[];summary=[];fixed_rows=[];verification=[]
    contexts=[('GSE166466','invasion'),('GSE149763','invasion'),('GSE147009','invasion'),('GSE166466','late_minus_early')]
    for label,context in contexts:
        c,e,s=cohorts[label]
        cc,X,contrast=primary(c,s) if context=='invasion' else stage(c,s)
        ranks=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True)
        fixed=-ranks.loc[targets].mean(axis=0)
        assert np.allclose(fixed,s.loc[cc.index,'Joint_program'])
        fit=fit_many(fixed.to_numpy()[:,None],X,contrast).iloc[0]
        independent=sm.OLS((fixed-fixed.mean())/fixed.std(ddof=1),X).fit(cov_type='HC3',use_t=True).t_test(contrast)
        assert np.isclose(fit.beta,float(np.asarray(independent.effect).item()))
        assert np.isclose(fit.p,float(np.asarray(independent.pvalue).item()),rtol=1e-9,atol=1e-12)
        verification.append({'context':label+'_'+context,'independent_statsmodels_match':True})
        fixed_rows.append({'cohort':label,'context':context,'kind':'Frozen_10_gene_program','n':len(cc),'measured':10,**fit.to_dict()})
        broad_have=[g for g in pool if g in ranks.index]
        broad=-ranks.loc[broad_have].mean(axis=0)
        broad_fit=fit_many(broad.to_numpy()[:,None],X,contrast).iloc[0]
        fixed_rows.append({'cohort':label,'context':context,'kind':'Nonoverlapping_broad_negative_score','n':len(cc),'measured':len(broad_have),**broad_fit.to_dict()})
        for family,sets in panels.items():
            ys=[];coverage=[]
            for genes in sets:
                have=[g for g in genes if g in ranks.index]
                coverage.append(len(have))
                ys.append(-ranks.loc[have].mean(axis=0).to_numpy() if len(have)>=8 else np.full(len(cc),np.nan))
            fits=fit_many(np.asarray(ys).T,X,contrast);fits['draw']=np.arange(B);fits['cohort']=label;fits['context']=context;fits['family']=family;fits['measured_genes']=coverage
            panel_fits.append(fits)
            good=fits.loc[fits.valid];n=len(good)
            exceed=int(good.t.ge(fit.t).sum())
            summary.append({'cohort':label,'context':context,'family':family,'n':len(cc),'draws_attempted':B,'scorable_draws':n,'fixed_beta':fit.beta,'fixed_t':fit.t,'fixed_p':fit.p,'reference_median_t':good.t.median(),'reference_positive_nominal_fraction':((good.beta>0)&(good.p<.05)).mean(),'reference_t_below_fixed_fraction':good.t.lt(fit.t).mean(),'add_one_reference_tail':(1+exceed)/(n+1),'reference_t_q025':good.t.quantile(.025),'reference_t_q975':good.t.quantile(.975),'limitation':'conditional gene-set reference; not phenotype permutation; correlation not matched'})
        print('Completed',label,context,'frozen t',fit.t,flush=True)
    out=pd.DataFrame(summary);out['holm_12_descriptive_reference_tails']=multipletests(out.add_one_reference_tail,method='holm')[1]
    out.to_csv(R/'program_distinctiveness_reference_summary.tsv',sep='\t',index=False)
    pd.concat(panel_fits,ignore_index=True).to_csv(R/'program_reference_all_fits.tsv.gz',sep='\t',index=False,compression='gzip')
    fixed_out=pd.DataFrame(fixed_rows);fixed_out.to_csv(R/'program_vs_broad_negative_score.tsv',sep='\t',index=False)
    # Estimate uncertainty without recoding a lack of significance as absence.
    tissue=pd.read_csv(R/'joint_program_all_tissue_cohorts.tsv',sep='\t')
    reference=float(tissue.loc[tissue.cohort.eq('GSE166466'),'beta'].iloc[0])
    rna=tissue[tissue.cohort.isin(['GSE149763','GSE147009','GSE223804'])].copy()
    rna['interval_contains_zero']=(rna.ci_low<=0)&(rna.ci_high>=0)
    rna['interval_contains_Clariom_point']=(rna.ci_low<=reference)&(rna.ci_high>=reference)
    rna['Clariom_point_reference_only']=reference
    rna['interpretation']='unconfirmed; interval overlap is not equivalence or a between-study test'
    rna.to_csv(R/'RNA_program_compatible_effect_ranges.tsv',sep='\t',index=False)
    pooled=pd.DataFrame([random_effects(rna[rna.cohort.isin(['GSE149763','GSE147009'])],'strictly_new_RNA_only'),random_effects(rna,'all_RNA_includes_prior_posthoc_context')])
    pooled.to_csv(R/'RNA_program_exploratory_pooled_estimation.tsv',sep='\t',index=False)
    checks={'seed':SEED,'draws_per_family':B,'families':3,'contexts':4,'all_pool':len(allpool),'negative_qualified_including_targets':len(negatives),'negative_reference_pool':len(pool),'program_targets':targets,'all_panels_unique_within_and_exclude_targets':True,'independent_statistical_checks':verification,'plan_sha256':hashlib.sha256((ROOT/'PROGRAM_DISTINCTIVENESS_v1.md').read_bytes()).hexdigest(),'meta_input_sha256':hashlib.sha256((R/'meta_discovery_all_genes.tsv').read_bytes()).hexdigest(),'completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'posthoc_after_original_outcomes_known':True}
    (R/'program_distinctiveness_run_record.json').write_text(json.dumps(checks,indent=2))
    fig,axs=plt.subplots(2,2,figsize=(11.8,7.5),layout='constrained')
    colors=['#7893a4','#c3a66b','#a56955']
    for ax,(label,context) in zip(axs.flat,contexts):
        for family,color in zip(panels,colors):
            row=out[(out.cohort==label)&(out.context==context)&(out.family==family)].iloc[0]
            data=next(d for d in panel_fits if d.cohort.iloc[0]==label and d.context.iloc[0]==context and d.family.iloc[0]==family)
            ax.hist(data.loc[data.valid,'t'],bins=35,density=True,histtype='step',lw=1.6,color=color,label=family.replace('_',' '))
        ax.axvline(row.fixed_t,color='#233e50',lw=1.8,ls='--',label='Frozen 10-gene score')
        ax.set(title=label+' | '+context.replace('_',' '),xlabel='Signed HC3 t statistic',ylabel='Reference density')
        if ax is axs.flat[0]:ax.legend(fontsize=7,frameon=False)
    fig.savefig(F/'16_program_reference_comparisons.png',dpi=300);fig.savefig(F/'16_program_reference_comparisons.pdf');plt.close(fig)
    print(out.to_string(index=False));print('\nRNA exploratory pooling:\n',pooled.to_string(index=False))

if __name__=='__main__':main()

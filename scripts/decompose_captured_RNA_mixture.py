"""Descriptive exact decomposition in captured single-cell RNA, not causal bulk."""
import json,hashlib,itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import ROOT,P,R,F
from extension_utils import common_universe,frozen_targets

def main():
    summary=pd.read_csv(R/'GSE208653_specimen_compartment_programs.tsv',sep='\t')
    counts=pd.read_csv(P/'GSE208653_compartment_pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
    donors=summary.loc[summary.stage.isin(['HSIL','SCC']),'donor'].drop_duplicates().tolist()
    stage=summary.drop_duplicates('donor').set_index('donor').stage
    sizes=summary[summary.donor.isin(donors)&~summary.compartment.eq('All cells')].pivot(index='compartment',columns='donor',values='cells')
    eligible=sizes.index[sizes.reindex(columns=donors).ge(30).all(axis=1)].tolist()
    assert 'Epithelial' in eligible and len(eligible)>=2
    old,targets,direction=frozen_targets();universe=counts.index.intersection(common_universe())
    def score(profile):
        rank=pd.Series(profile,index=counts.index).loc[universe].rank(pct=True)
        return -float(rank.loc[targets].mean())
    profiles={};weights={};coverage=[];observed=[]
    for donor in donors:
        columns=[donor+'|'+x for x in eligible]
        n=counts[columns].to_numpy(float);sums=n.sum(axis=0);w=sums/sums.sum();e=n/sums*1e6
        reconstructed=e@w;direct=n.sum(axis=1)/n.sum()*1e6
        assert np.allclose(reconstructed,direct,rtol=1e-12,atol=1e-10)
        profiles[donor]=e;weights[donor]=w
        coverage.append({'donor':donor,'stage':stage[donor],'eligible_compartments':';'.join(eligible),'retained_cell_fraction':float(sizes.loc[eligible,donor].sum()/summary.loc[(summary.donor==donor)&summary.compartment.eq('All cells'),'cells'].iloc[0]),'retained_captured_UMI_fraction':float(n.sum()/counts[donor+'|All cells'].sum()),'reconstruction_max_absolute_CPM_error':float(abs(reconstructed-direct).max())})
        for comp,weight in zip(eligible,w):observed.append({'donor':donor,'stage':stage[donor],'compartment':comp,'captured_RNA_weight':weight})
    rows=[]
    for h,s in itertools.product([d for d in donors if stage[d]=='HSIL'],[d for d in donors if stage[d]=='SCC']):
        eh,es=profiles[h],profiles[s];wh,ws=weights[h],weights[s]
        hh,ss,sh,hs=score(eh@wh),score(es@ws),score(eh@ws),score(es@wh)
        composition=((sh-hh)+(ss-hs))/2
        profile=((hs-hh)+(ss-sh))/2
        total=ss-hh
        assert np.isclose(composition+profile,total,atol=1e-14)
        rows.append({'HSIL_specimen':h,'SCC_specimen':s,'score_HSIL_observed_restricted':hh,'score_SCC_observed_restricted':ss,'observed_restricted_difference':total,'captured_RNA_composition_component':composition,'within_compartment_profile_component':profile,'exact_identity_error':composition+profile-total,'interpretation':'crossing pairs share specimens; descriptive only; not a fraction of actual bulk change'})
    common=np.mean(list(weights.values()),axis=0)
    standardized=[]
    for d in donors:standardized.append({'donor':d,'stage':stage[d],'score_observed_restricted':score(profiles[d]@weights[d]),'score_at_common_captured_RNA_weights':score(profiles[d]@common)})
    out=pd.DataFrame(rows);out.to_csv(R/'captured_RNA_mixture_pair_decomposition.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(R/'captured_RNA_mixture_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(observed).to_csv(R/'captured_RNA_mixture_weights.tsv',sep='\t',index=False)
    std=pd.DataFrame(standardized);std.to_csv(R/'captured_RNA_mixture_standardized_scores.tsv',sep='\t',index=False)
    (R/'captured_RNA_mixture_record.json').write_text(json.dumps({'specimens':4,'independent_HSIL_specimens':2,'independent_SCC_specimens':2,'cross_pairs_are_not_independent':True,'eligible_compartments':eligible,'weight_type':'fraction of captured UMI counts, not true tissue cell fractions','all_exact_decomposition_and_reconstruction_checks_passed':True,'causal_interpretation':False,'plan_sha256':hashlib.sha256((ROOT/'SINGLE_CELL_MIXTURE_v1.md').read_bytes()).hexdigest()},indent=2))
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    x=np.arange(len(out));width=.24
    for i,(col,label,color) in enumerate([('observed_restricted_difference','Observed restricted mixture','#3f6478'),('captured_RNA_composition_component','Captured-RNA weights','#bf9a60'),('within_compartment_profile_component','Within-compartment profile','#a66051')]):
        axs[0].bar(x+(i-1)*width,out[col],width,label=label,color=color)
    axs[0].axhline(0,color='grey',lw=.8);axs[0].set(xticks=x,xticklabels=['684→686','684→687','685→686','685→687'],ylabel='Difference in fixed signed rank score',title='A  Synthetic-mixture decomposition')
    axs[0].legend(frameon=False,fontsize=7)
    for _,row in std.iterrows():
        axs[1].plot([0,1],[row.score_observed_restricted,row.score_at_common_captured_RNA_weights],'-o',color='#c4a15d' if row.stage=='HSIL' else '#a76050')
        dy=8 if row.donor.endswith('684') else -8 if row.donor.endswith('685') else 0
        axs[1].annotate(row.donor[-3:],(1,row.score_at_common_captured_RNA_weights),xytext=(6,dy),textcoords='offset points',fontsize=8)
    axs[1].set(xticks=[0,1],xticklabels=['Observed weights','Common weights'],ylabel='Fixed signed rank score',title='B  One point per specimen');axs[1].set_xlim(-.12,1.28)
    fig.savefig(F/'17_captured_RNA_mixture_sensitivity.png',dpi=300);fig.savefig(F/'17_captured_RNA_mixture_sensitivity.pdf');plt.close(fig)
    print('Eligible compartments',eligible);print(pd.DataFrame(coverage).to_string(index=False));print(out.to_string(index=False));print(std.to_string(index=False))

if __name__=='__main__':main()

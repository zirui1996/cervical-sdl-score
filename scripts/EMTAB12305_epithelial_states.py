"""Reuse the old independent marker-state rule in the new primary specimens."""
import gc,gzip,hashlib,json
import numpy as np
import pandas as pd
from scipy.io import mmread
from run_analysis import ROOT,RAW,P,R
from epithelial_context_checks import SETS
from extension_utils import extra_scores

def main():
    cm=pd.read_csv(P/'EMTAB12305_cell_QC_marker_compartments.tsv.gz',sep='\t')
    folder=RAW/'E-MTAB-12305';pb=[];rows=[];reviews=[]
    for donor in ['H1','H2','T1','T2','T3','T4']:
        cells=cm[cm.specimen.eq(donor)&cm.marker_compartment.eq('Epithelial')];stage=cells.stage.iloc[0]
        feat=pd.read_csv(folder/(donor+'features.tsv.gz'),sep='\t',header=None);symbols=feat.iloc[:,1].astype(str).to_numpy()
        bc=pd.read_csv(folder/(donor+'barcodes.tsv.gz'),sep='\t',header=None).iloc[:,0].astype(str)
        ix=pd.Index(bc).get_indexer(cells.barcode);assert (ix>=0).all()
        with gzip.open(folder/(donor+'matrix.mtx.gz'),'rb') as f:X=mmread(f).tocsr().astype(np.int32)
        X=X[:,ix].tocsr();lib=np.asarray(X.sum(axis=0)).ravel();scores=[];det=[]
        for name,markers in SETS.items():
            vals=[np.asarray(X[symbols==g].sum(axis=0)).ravel() for g in markers]
            scores.append(np.mean([np.log1p(v/lib*1e4) for v in vals],axis=0));det.append(np.sum([v>0 for v in vals],axis=0))
        score=np.asarray(scores);det=np.asarray(det);order=np.argsort(score,axis=0);best=order[-1];second=order[-2];k=np.arange(len(ix))
        good=(score[best,k]>=.3)&(det[best,k]>=2)&(score[best,k]-score[second,k]>=.15)
        names=np.asarray(list(SETS),dtype=object)[best];names[~good]='Mixed_unspecified'
        for state in sorted(set(names)):
            mask=names==state;counts=np.asarray(X[:,mask].sum(axis=1)).ravel()
            key=donor+'|'+state;vec=pd.DataFrame({'gene':symbols,'count':counts}).groupby('gene')['count'].sum();pb.append(vec.rename(key))
            rows.append({'key':key,'specimen':donor,'stage':stage,'marker_state':state,'cells':int(mask.sum()),'fraction_of_epithelial':float(mask.mean())})
            for j,name in enumerate(SETS):reviews.append({'specimen':donor,'stage':stage,'state':state,'marker_set':name,'mean_marker_score':float(score[j,mask].mean())})
        print(donor,pd.Series(names).value_counts().to_dict(),flush=True);del X;gc.collect()
    counts=pd.concat(pb,axis=1).fillna(0).astype(np.int64)
    counts.to_csv(P/'EMTAB12305_epithelial_state_pseudobulk_counts.tsv.gz',sep='\t',compression='gzip')
    e=np.log2(counts.div(counts.sum(),axis=1)*1e6+1);s,cov=extra_scores(e)
    out=pd.DataFrame(rows).set_index('key').join(s[['Joint_program','Composite']]);out.loc[out.cells<30,'Joint_program']=np.nan
    out.reset_index().to_csv(R/'EMTAB12305_epithelial_marker_states.tsv',sep='\t',index=False)
    pd.DataFrame(reviews).to_csv(R/'EMTAB12305_epithelial_state_marker_review.tsv',sep='\t',index=False)
    (R/'EMTAB12305_marker_states_record.json').write_text(json.dumps({'uses_exact_prior_marker_sets_and_thresholds':True,'posthoc_after_new_program_results':True,'not_histological_or_malignant_identity':True,'plan_sha256':hashlib.sha256((ROOT/'EMTAB12305_MARKER_STATES_v1.md').read_bytes()).hexdigest()},indent=2))
    print(out.to_string(),flush=True)

if __name__=='__main__':main()

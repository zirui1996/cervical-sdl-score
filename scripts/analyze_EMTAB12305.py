"""Official count-level localization extension; specimen-level description only."""
import gc,gzip,hashlib,json
import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import ROOT,RAW,P,R,F
from extension_utils import extra_scores,frozen_targets
from analyze_singlecell_extension import MARKERS

def main():
    folder=RAW/'E-MTAB-12305'
    sdrf=pd.read_csv(folder/'E-MTAB-12305.sdrf.txt',sep='\t')
    sdrf['specimen']=sdrf['Derived Array Data File'].str.replace('barcodes.tsv.gz','',regex=False)
    sdrf['stage']=sdrf['Characteristics[sampling site]'].map({'high squamous intraepithelial lesion':'HSIL','neoplasm':'SCC','normal tissue adjacent to neoplasm':'Adjacent_normal','metastatic lymph node':'Metastatic_lymph_node'})
    assert sdrf.stage.notna().all()
    metadata=sdrf[['specimen','Source Name','Characteristics[age]','Characteristics[disease]','Characteristics[sampling site]','stage']]
    metadata.to_csv(P/'EMTAB12305_specimen_metadata.tsv',sep='\t',index=False)
    chosen=metadata[~metadata.stage.eq('Metastatic_lymph_node')].copy()
    ordering={x:i for i,x in enumerate(['H1','H2','T1','T2','T3','T4','N1','N2','N3'])}
    chosen=chosen.assign(order=chosen.specimen.map(ordering)).sort_values('order')
    assert chosen.stage.value_counts().to_dict()=={'SCC':4,'Adjacent_normal':3,'HSIL':2}
    old,meta,directions=frozen_targets();targets=list(directions)
    review=sorted(set(sum(MARKERS.values(),[])))
    qc=[];records=[];localization=[];marker_review=[];pseudobulk=[];annotations=[]
    for _,sample in chosen.iterrows():
        donor=sample.specimen;stage=sample.stage
        paths={kind:folder/(donor+suffix) for kind,suffix in [('features','features.tsv.gz'),('barcodes','barcodes.tsv.gz'),('matrix','matrix.mtx.gz')]}
        for p in paths.values():assert p.exists(),str(p)
        features=pd.read_csv(paths['features'],sep='\t',header=None)
        symbols=features.iloc[:,1].astype(str).to_numpy();barcodes=pd.read_csv(paths['barcodes'],sep='\t',header=None).iloc[:,0].astype(str).to_numpy()
        with gzip.open(paths['matrix'],'rb') as f:X=mmread(f).tocsr()
        assert X.shape==(len(features),len(barcodes))
        assert X.data.min()>=0 and np.allclose(X.data,np.round(X.data)),'Count units not verified'
        X=X.astype(np.int32)
        total=np.asarray(X.sum(axis=0)).ravel();detected=np.asarray((X>0).sum(axis=0)).ravel()
        mt=np.asarray(X[np.char.startswith(symbols.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(total,1)
        keep=(detected>=200)&(detected<=7500)&(mt<=.2)&(total>0)
        qc.append({'specimen':donor,'stage':stage,'input_cells':X.shape[1],'retained_cells':int(keep.sum()),'excluded_cells':int((~keep).sum()),'nonnegative_integer_counts_verified':True})
        X=X[:,keep].tocsr();barcodes=barcodes[keep];lib=total[keep]
        index={g:np.flatnonzero(symbols==g) for g in set(review+targets)}
        raw={g:np.asarray(X[ix].sum(axis=0)).ravel() if len(ix) else np.zeros(X.shape[1]) for g,ix in index.items()}
        norm={g:np.log1p(v/np.maximum(lib,1)*1e4) for g,v in raw.items()}
        names=list(MARKERS);score=np.zeros((len(names),X.shape[1]));present=np.zeros_like(score)
        for j,lineage in enumerate(names):
            have=[g for g in MARKERS[lineage] if len(index[g])]
            score[j]=np.mean([norm[g] for g in have],axis=0);present[j]=np.sum([raw[g]>0 for g in have],axis=0)
        order=np.argsort(score,axis=0);first=order[-1];second=order[-2];ix=np.arange(X.shape[1])
        conf=(score[first,ix]>=.3)&(score[first,ix]-score[second,ix]>=.15)&(present[first,ix]>=2)
        labels=np.asarray(names,dtype=object)[first];labels[~conf]='Unassigned ambiguous'
        annotations.append(pd.DataFrame({'specimen':donor,'stage':stage,'barcode':barcodes,'marker_compartment':labels,'total_counts':lib,'detected_genes':detected[keep],'mito_fraction':mt[keep]}))
        for comp in ['All cells']+sorted(set(labels)):
            mask=np.ones(len(labels),bool) if comp=='All cells' else labels==comp
            sums=np.asarray(X[:,mask].sum(axis=1)).ravel()
            pb=pd.DataFrame({'gene':symbols,'count':sums}).groupby('gene')['count'].sum()
            key=donor+'|'+comp;pseudobulk.append(pb.rename(key))
            records.append({'specimen':donor,'stage':stage,'compartment':comp,'cells':int(mask.sum()),'fraction_of_QC_cells':float(mask.mean()),'pseudobulk_total_counts':int(sums.sum()),'key':key})
            for gene in targets:localization.append({'specimen':donor,'stage':stage,'compartment':comp,'gene':gene,'cells':int(mask.sum()),'detected_fraction':float((raw[gene][mask]>0).mean()),'mean_log1p_CP10k':float(norm[gene][mask].mean()),'pseudobulk_counts':int(raw[gene][mask].sum())})
            if comp!='All cells':
                for gene in review:marker_review.append({'specimen':donor,'stage':stage,'compartment':comp,'gene':gene,'cells':int(mask.sum()),'mean_log1p_CP10k':float(norm[gene][mask].mean()),'detected_fraction':float((raw[gene][mask]>0).mean())})
        print(donor,stage,qc[-1]['input_cells'],'->',qc[-1]['retained_cells'],pd.Series(labels).value_counts().to_dict(),flush=True)
        del X,features,raw,norm;gc.collect()
    counts=pd.concat(pseudobulk,axis=1).fillna(0).astype(np.int64)
    counts.to_csv(P/'EMTAB12305_compartment_pseudobulk_counts.tsv.gz',sep='\t',compression='gzip')
    pd.concat(annotations,ignore_index=True).to_csv(P/'EMTAB12305_cell_QC_marker_compartments.tsv.gz',sep='\t',index=False,compression='gzip')
    e=np.log2(counts.div(counts.sum(axis=0),axis=1)*1e6+1);s,cov=extra_scores(e)
    out=pd.DataFrame(records).set_index('key').join(s[['Joint_program','Composite']]);out.loc[out.cells<30,'Joint_program']=np.nan
    out.reset_index().to_csv(R/'EMTAB12305_specimen_compartment_programs.tsv',sep='\t',index=False)
    pd.DataFrame(qc).to_csv(R/'EMTAB12305_QC_summary.tsv',sep='\t',index=False)
    pd.DataFrame(localization).to_csv(R/'EMTAB12305_gene_localization.tsv',sep='\t',index=False)
    pd.DataFrame(marker_review).to_csv(R/'EMTAB12305_marker_review.tsv',sep='\t',index=False)
    cov.to_csv(R/'EMTAB12305_program_coverage.tsv',sep='\t',index=False)
    epi=out[out.compartment.eq('Epithelial')]
    valid=epi[epi.cells.ge(30)&epi.Joint_program.notna()]
    supported=set(valid.loc[valid.stage.isin(['HSIL','SCC']),'specimen'])==set(chosen.loc[chosen.stage.isin(['HSIL','SCC']),'specimen'])
    (R/'EMTAB12305_analysis_summary.json').write_text(json.dumps({'specimens_analyzed':len(chosen),'adjacent_normal_specimens':3,'HSIL_specimens':2,'SCC_specimens':4,'lymph_node_excluded':True,'input_cells':sum(x['input_cells'] for x in qc),'retained_cells':sum(x['retained_cells'] for x in qc),'all_primary_specimens_have_30_epithelial_cells_and_scorable_program':supported,'patient_pairing_not_inferred_from_equal_ages':True,'not_malignant_cell_identity':True,'no_cellwise_inference':True,'source':'EBI-deposited sparse integer-count matrices; see each transfer provenance, including any hash-verified mirror; not author log-normalized dense CSV','plan_sha256':hashlib.sha256((ROOT/'REFERENCE_METHOD_TRANSFER_v1.md').read_bytes()).hexdigest()},indent=2))
    fig,axs=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
    colors={'Adjacent_normal':'#74939b','HSIL':'#c8a04e','SCC':'#a95e51'}
    for ax,comp in zip(axs,['All cells','Epithelial']):
        for j,(group,color) in enumerate(colors.items()):
            part=out[out.compartment.eq(comp)&out.stage.eq(group)&out.cells.ge(30)]
            x=np.linspace(j-.13,j+.13,len(part))
            ax.scatter(x,part.Joint_program,color=color,s=42)
            for xx,(_,row) in zip(x,part.iterrows()):ax.annotate(row.name.split('|')[0],(xx,row.Joint_program),xytext=(3,4),textcoords='offset points',fontsize=7)
        ax.set(xticks=range(3),xticklabels=['Adjacent normal','HSIL','SCC'],ylabel='Frozen signed mean fractional rank',title=comp+' | one point per specimen',xlim=(-.35,2.45))
    fig.savefig(F/'18_EMTAB12305_fixed_program.png',dpi=300);fig.savefig(F/'18_EMTAB12305_fixed_program.pdf');plt.close(fig)
    print('\nSpecimen epithelial score:\n',epi[['stage','cells','Joint_program','Composite']].to_string(),flush=True)

if __name__=='__main__':main()

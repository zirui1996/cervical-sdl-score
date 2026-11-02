"""Memory-bounded specimen-by-compartment localization; no cell-wise inference."""
from pathlib import Path
import gc,gzip,json,re
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import RAW,P,R,F,dump
from extension_utils import frozen_targets,extra_scores

MARKERS={
 'Epithelial':['EPCAM','TACSTD2','KRT5','KRT14','KRT17','KRT19','KRT8','KRT18','KRT7','KRT13','KRT4'],
 'Immune':['PTPRC','CD3D','CD3E','TRAC','CD79A','MS4A1','LST1','TYROBP','FCER1G','NKG7','GNLY'],
 'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA'],
 'Endothelial':['PECAM1','VWF','EMCN','KDR','PLVAP'],
 'Smooth muscle pericyte':['ACTA2','TAGLN','MYH11','RGS5','CSPG4'],
 'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],
 'Erythroid':['HBB','HBA1','HBA2','ALAS2','GYPA']}

def main():
 old,meta,directions=frozen_targets();targets=list(directions);review=sorted(set(sum(MARKERS.values(),[]))|{'FCGR3B','CSF3R','S100A8','S100A9'})
 qc=[];summaries=[];genes=[];markerrows=[];pseudobulk=[];cellmeta=[]
 files=sorted(RAW.glob('GSM636068*.matrix.mtx.gz'));assert len(files)==9
 for f in files:
  donor=f.name.split('_')[0];base=f.name[:-len('.matrix.mtx.gz')]
  stage='ADC' if '_ADC_' in base else 'SCC' if '_SCC_' in base else 'HSIL' if '_HSIL_' in base else 'HPV_negative_normal' if '_HPV_NEG_' in base else 'HPV_positive_normal'
  features=pd.read_csv(RAW/(base+'.features.tsv.gz'),sep='\t',header=None);symbols=features.iloc[:,1].astype(str).to_numpy()
  barcodes=pd.read_csv(RAW/(base+'.barcodes.tsv.gz'),sep='\t',header=None).iloc[:,0].astype(str).to_numpy()
  with gzip.open(f,'rb') as h:X=mmread(h).tocsr().astype(np.int32)
  assert X.shape==(len(features),len(barcodes)) and X.data.min()>=0
  total=np.asarray(X.sum(axis=0)).ravel();detected=np.asarray((X>0).sum(axis=0)).ravel();mt=np.asarray(X[np.char.startswith(symbols.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(total,1)
  keep=(detected>=200)&(detected<=7500)&(mt<=.20)&(total>0)
  qc.append({'donor':donor,'stage':stage,'input_cells':X.shape[1],'retained_cells':int(keep.sum()),'excluded_cells':int((~keep).sum()),'median_detected_genes_retained':float(np.median(detected[keep])),'median_mito_fraction_retained':float(np.median(mt[keep]))})
  X=X[:,keep].tocsr();barcodes=barcodes[keep];lib=total[keep];geneindex={g:np.flatnonzero(symbols==g) for g in set(review+targets)}
  selected={}
  for g,ix in geneindex.items():selected[g]=np.asarray(X[ix].sum(axis=0)).ravel() if len(ix) else np.zeros(X.shape[1])
  norm={g:np.log1p(v/lib*1e4) for g,v in selected.items()}
  names=list(MARKERS);score=np.zeros((len(names),X.shape[1]));present=np.zeros_like(score)
  for j,lineage in enumerate(names):
   have=[g for g in MARKERS[lineage] if len(geneindex[g])]
   score[j]=np.mean([norm[g] for g in have],axis=0);present[j]=np.sum([selected[g]>0 for g in have],axis=0)
  order=np.argsort(score,axis=0);first=order[-1];second=order[-2];ind=np.arange(X.shape[1]);conf=(score[first,ind]>=.3)&((score[first,ind]-score[second,ind])>=.15)&(present[first,ind]>=2)
  lineage=np.array(names,dtype=object)[first];lineage[~conf]='Unassigned ambiguous'
  # Save only necessary per-cell QC and broad marker label; raw matrices remain
  # separate, avoiding a huge in-memory concatenation or artificial replication.
  cm=pd.DataFrame({'donor':donor,'stage':stage,'barcode':barcodes,'marker_compartment':lineage,'total_counts':lib,'detected_genes':detected[keep],'mito_fraction':mt[keep]});cellmeta.append(cm)
  for group in ['All cells']+sorted(set(lineage)):
   ix=np.ones(len(lineage),bool) if group=='All cells' else lineage==group
   counts=np.asarray(X[:,ix].sum(axis=1)).ravel();pb=pd.DataFrame({'gene':symbols,'count':counts}).groupby('gene').sum()['count'];key=donor+'|'+group
   pseudobulk.append(pb.rename(key))
   summaries.append({'donor':donor,'stage':stage,'compartment':group,'cells':int(ix.sum()),'fraction_of_QC_cells':float(ix.mean()),'pseudobulk_total_counts':int(counts.sum()),'pseudobulk_key':key})
   for g in targets:genes.append({'donor':donor,'stage':stage,'compartment':group,'gene':g,'cells':int(ix.sum()),'detected_fraction':float((selected[g][ix]>0).mean()),'mean_log1p_CP10k':float(norm[g][ix].mean()),'pseudobulk_counts':int(selected[g][ix].sum())})
   if group!='All cells':
    for g in review:markerrows.append({'donor':donor,'stage':stage,'compartment':group,'gene':g,'cells':int(ix.sum()),'mean_log1p_CP10k':float(norm[g][ix].mean()),'detected_fraction':float((selected[g][ix]>0).mean())})
  print(donor,stage,qc[-1]['input_cells'],'->',qc[-1]['retained_cells'],pd.Series(lineage).value_counts().to_dict(),flush=True)
  del X,features,selected,norm;gc.collect()
 countmat=pd.concat(pseudobulk,axis=1).fillna(0).astype(np.int64);countmat.to_csv(P/'GSE208653_compartment_pseudobulk_counts.tsv.gz',sep='\t',compression='gzip')
 pd.concat(cellmeta,ignore_index=True).to_csv(P/'GSE208653_cell_QC_marker_compartments.tsv.gz',sep='\t',index=False,compression='gzip')
 e=np.log2(countmat.div(countmat.sum(axis=0),axis=1)*1e6+1);s,cov=extra_scores(e)
 summary=pd.DataFrame(summaries).set_index('pseudobulk_key');summary=summary.join(s[['Composite','Joint_program']]);summary.loc[summary.cells<30,'Joint_program']=np.nan
 summary.reset_index().to_csv(R/'GSE208653_specimen_compartment_programs.tsv',sep='\t',index=False)
 pd.DataFrame(genes).to_csv(R/'GSE208653_specimen_gene_localization.tsv',sep='\t',index=False);pd.DataFrame(markerrows).to_csv(R/'GSE208653_compartment_marker_review.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(R/'GSE208653_QC_summary.tsv',sep='\t',index=False);cov.to_csv(R/'GSE208653_pseudobulk_score_coverage.tsv',sep='\t',index=False)
 dump('GSE208653_extension_summary.json',{'specimens':len(qc),'input_cells':int(sum(x['input_cells'] for x in qc)),'retained_cells':int(sum(x['retained_cells'] for x in qc)),'QC':'200-7500 genes; mitochondrial fraction <=0.20; positive counts','classification':'fixed canonical-marker broad compartments with ambiguity retained; not malignant-cell calls','inference':'specimen-level descriptive only; 2 HSIL versus 2 SCC; ADC excluded from squamous comparison','targets':targets})
 # Donor-level points, never cell-wise significance stars.
 colors={'HPV_negative_normal':'#6c8796','HPV_positive_normal':'#7ca77c','HSIL':'#d4a350','SCC':'#a95d50','ADC':'#83729b'}
 fig,axs=plt.subplots(1,2,figsize=(11,4.1))
 for ax,comp in zip(axs,['All cells','Epithelial']):
  v=summary[summary.compartment.eq(comp)&summary.cells.ge(30)]
  for k,stage in enumerate(colors):
   vals=v[v.stage.eq(stage)].Joint_program;xx=np.linspace(k-.09,k+.09,len(vals)) if len(vals)>1 else [k]
   ax.scatter(xx,vals,color=colors[stage],s=45)
  ax.set(xticks=range(5),xticklabels=['HPV- normal','HPV+ normal','HSIL','SCC','ADC'],ylabel='Fixed program signed mean fractional rank',title=comp+' | one point per specimen');ax.tick_params(axis='x',rotation=30)
 fig.tight_layout();fig.savefig(F/'12_singlecell_specimen_program.png',dpi=300);fig.savefig(F/'12_singlecell_specimen_program.pdf');plt.close(fig)

if __name__=='__main__':main()

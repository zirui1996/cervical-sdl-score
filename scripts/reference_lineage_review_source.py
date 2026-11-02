"""Manual broad-lineage review based on canonical markers, independent of candidate genes."""
from pathlib import Path
import json,os
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data/processed';R=ROOT/'results';F=ROOT/'figures'
a=ad.read_h5ad(P/'GSE171894_analyzed.h5ad')
review={
 '0':('Epithelial','SPRR3 TACSTD2 KRT13 CLDN7; low PTPRC/FCGR3B/CSF3R'),
 '1':('Epithelial','SERPINB3/4 DUOX2 KRT6A KRT13 TACSTD2; S100A8/9 alone are not neutrophil-specific'),
 '2':('Epithelial','KRT5 KRT13 SFN PERP TACSTD2'),
 '3':('Erythroid / mixed','HBB HBA1 HBA2 enrichment with residual epithelial/immune markers; mixed identity retained'),
 '4':('Myeloid','LYZ HLA-DRA CD74 FCER1G IL1B'),
 '5':('T cells','TRAC CD3D CD3E TRBC2 ICOS'),
 '6':('T cells','TRAC CD3D CD8A CCL5; cytotoxic T program'),
 '7':('T cells','TRAC CD3D IL7R'),
 '8':('Plasma cells','MZB1 JCHAIN immunoglobulins CD79A'),
 '9':('Epithelial','KRT5 KRT13 SFN SERPINB3; low PTPRC'),
 '10':('B cells','MS4A1 CD79A CD74 HLA-DRA'),
 '11':('Dendritic (pDC-like)','LILRA4 CLEC4C TCF4 GZMB; low MS4A1/CD79A and MZB1'),
 '12':('T/NK cytotoxic','GNLY NKG7 CCL5 with appreciable TRAC/CD3D; broad mixed cytotoxic lineage retained'),
 '13':('NK cells','GNLY KLRD1 FCER1G IL2RB; low CD3D'),
 '14':('Epithelial','SPRR2A/D SERPINB4 PI3 SLPI TACSTD2; absent canonical granulocyte markers')}
assert set(a.obs.leiden.astype(str).unique())==set(review),'Cluster identities changed; review is required again.'
old=a.obs.lineage.astype(str).copy();a.obs['lineage']=a.obs.leiden.astype(str).map({k:v[0] for k,v in review.items()}).astype('category')
pd.DataFrame([{'cluster':k,'reviewed_lineage':v[0],'marker_basis':v[1],'cells':int(a.obs.leiden.astype(str).eq(k).sum())} for k,v in review.items()]).to_csv(R/'singlecell_manual_annotation_review.tsv',sep='\t',index=False)
a.uns['manual_lineage_review']='Canonical-marker review after automated broad annotation; target candidate expression was not used as an annotation feature. Clusters 0/1/14 corrected from neutrophil to epithelial; cluster 11 pDC-like; cluster 3 erythroid/mixed; cluster 12 T/NK.'
tmp=P/'GSE171894_reviewed.complete.h5ad';a.write_h5ad(tmp,compression='gzip');check=ad.read_h5ad(tmp);assert check.shape==a.shape;del check;os.replace(tmp,P/'GSE171894_analyzed.h5ad')
a.obs.assign(UMAP1=a.obsm['X_umap'][:,0],UMAP2=a.obsm['X_umap'][:,1]).to_csv(P/'singlecell_cell_annotations.tsv.gz',sep='\t',compression='gzip')
counts=pd.crosstab(a.obs.donor,a.obs.lineage);counts.to_csv(R/'singlecell_lineage_counts.tsv',sep='\t')
genes=[r['gene'] for r in json.loads((R/'candidate_lock.json').read_text())['candidates']];rows=[]
for (donor,lineage),ix in a.obs.groupby(['donor','lineage'],observed=True).groups.items():
 for gene in genes:
  v=a[ix,gene].X.toarray().ravel();rows.append({'donor':donor,'lineage':str(lineage),'gene':gene,'n_cells':len(v),'mean_author_normalized_expression':float(v.mean()),'detection_fraction':float((v>0).mean())})
pd.DataFrame(rows).to_csv(R/'singlecell_candidate_by_donor_lineage.tsv',sep='\t',index=False)
summary=json.loads((R/'singlecell_summary.json').read_text());summary['annotation']='canonical-marker-reviewed broad lineages; malignant identity unconfirmed; mixed lineages retained';summary['manual_reviewed']=True;(R/'singlecell_summary.json').write_text(json.dumps(summary,indent=2))
fig,axs=plt.subplots(1,3,figsize=(16,5));u=a.obsm['X_umap']
for ax,col in zip(axs[:2],['lineage','donor']):
 for i,category in enumerate(a.obs[col].astype('category').cat.categories):
  mask=a.obs[col].eq(category).values;ax.scatter(u[mask,0],u[mask,1],s=2,color=plt.get_cmap('tab20')(i),label=str(category),rasterized=True)
 ax.set(title=col,xlabel='UMAP1',ylabel='UMAP2');ax.legend(fontsize=6,markerscale=4,loc='best')
v=a[:,genes[0]].X.toarray().ravel();ix=np.argsort(v);im=axs[2].scatter(u[ix,0],u[ix,1],c=v[ix],s=2,cmap='magma',vmax=np.quantile(v,.99),rasterized=True);fig.colorbar(im,ax=axs[2],label='Author-normalized expression');axs[2].set(title=genes[0],xlabel='UMAP1',ylabel='UMAP2')
fig.suptitle('GSE171894: marker-reviewed descriptive cell lineages (4 donors)');fig.tight_layout();fig.savefig(F/'03_singlecell_overview.png',dpi=150);fig.savefig(F/'03_singlecell_overview.pdf');plt.close(fig)
print('Reviewed annotations:',counts.sum().to_dict());print('Cells with revised lineage:',int((old.values!=a.obs.lineage.astype(str).values).sum()))

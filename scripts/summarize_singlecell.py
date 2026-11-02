"""Donor-level descriptive localization using the previously analyzed public reference."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import ROOT,P,R,F,dump

initial=[r['gene'] for r in json.loads((R/'candidate_lock.json').read_text())['candidates']]
extended=json.loads((R/'RNAseq_candidate_lock.json').read_text())['genes'];genes=initial+extended
reference=P/'GSE171894_analyzed_reference.h5ad'
if not reference.exists():
 source=ROOT.parent/'cervical_ln_bioinformatics/data/processed/GSE171894_analyzed.h5ad'
 assert source.exists(),'Analyzed public single-cell reference must be supplied.'
 shutil.copyfile(source,reference)
a=ad.read_h5ad(reference);rows=[]
for (donor,lineage),ix in a.obs.groupby(['donor','lineage'],observed=True).groups.items():
 for gene in genes:
  if gene not in a.var_names:continue
  v=a[ix,gene].X.toarray().ravel()
  rows.append({'gene':gene,'donor':donor,'lineage':str(lineage),'cells':len(v),'mean_author_normalized_expression':float(v.mean()),'detection_fraction':float((v>0).mean())})
d=pd.DataFrame(rows);d.to_csv(R/'singlecell_donor_lineage_expression.tsv',sep='\t',index=False)
ag=d.groupby(['gene','lineage']).agg(donors=('donor','nunique'),cells=('cells','sum'),mean_of_donor_means=('mean_author_normalized_expression','mean'),mean_detection_fraction=('detection_fraction','mean')).reset_index();ag.to_csv(R/'singlecell_summary_by_lineage.tsv',sep='\t',index=False)
present=[g for g in genes if g in a.var_names];small=a[:,present].copy();small.write_h5ad(P/'GSE171894_selected_gene_localization.h5ad',compression='gzip')
record={'reference_dataset':'GSE171894','source_url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171894/suppl/GSE171894_RAW.tar','source_raw_sha256':hashlib.sha256((ROOT.parent/'cervical_ln_bioinformatics/data/raw/GSE171894_RAW.tar').read_bytes()).hexdigest() if (ROOT.parent/'cervical_ln_bioinformatics/data/raw/GSE171894_RAW.tar').exists() else None,'reference_sha256':hashlib.sha256(reference.read_bytes()).hexdigest(),'cells':a.n_obs,'donors':int(a.obs.donor.nunique()),'annotation':'reused provisional broad marker-based lineage labels; epithelial is not confirmed malignant','units':'author-normalized UMI values; precise normalization transform unspecified','absent_genes':sorted(set(genes)-set(a.var_names)),'inference':'descriptive donor summaries only; no cellular pseudoreplication tests'}
dump('singlecell_reference_provenance.json',record)
mat=ag.pivot(index='gene',columns='lineage',values='mean_of_donor_means').reindex(genes)
fig,ax=plt.subplots(figsize=(10,6));v=mat.values;scaled=v/(v.max(axis=1,keepdims=True)+1e-12);im=ax.imshow(scaled,aspect='auto',vmin=0,vmax=1,cmap='YlOrRd');ax.set(xticks=range(len(mat.columns)),xticklabels=mat.columns,yticks=range(len(mat.index)),yticklabels=mat.index,title='Descriptive lineage localization: donor-weighted mean, scaled within each gene');plt.setp(ax.get_xticklabels(),rotation=40,ha='right');fig.colorbar(im,ax=ax,label='Relative mean within gene');fig.tight_layout();fig.savefig(F/'06_singlecell_localization.png',dpi=150);fig.savefig(F/'06_singlecell_localization.pdf');plt.close(fig)
print(json.dumps(record,indent=2));print(ag[ag.gene.isin(extended)].to_string(index=False))

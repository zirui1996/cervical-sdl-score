import json
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from download_data import download,RAW
from run_analysis import R

name='c5.go.bp.v2025.1.Hs.symbols.gmt';url='https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/'+name
if not (RAW/name).exists():download((name,url))
d=pd.read_csv(R/'meta_discovery_all_genes.tsv',sep='\t');background=set(d.gene);qual=d[d.same_direction_all_three&d.combined_fdr.lt(.05)&d.mean_absolute_standardized_beta.ge(.5)]
program={x['gene'] for x in json.loads((R/'joint_program_lock.json').read_text())['candidates']}
foregrounds={'Frozen_ten_gene_program':program,'Qualified_positive':set(qual.loc[qual.direction.eq(1),'gene']),'Qualified_negative':set(qual.loc[qual.direction.eq(-1),'gene'])}
terms={}
for line in (RAW/name).read_text().splitlines():
 row=line.split('\t');genes=set(row[2:])&background
 if 10<=len(genes)<=500:terms[row[0]]=genes
rows=[]
for label,fg in foregrounds.items():
 for term,genes in terms.items():
  overlap=fg&genes;rows.append({'foreground':label,'term':term,'foreground_genes':len(fg),'term_background_genes':len(genes),'background_genes':len(background),'overlap':len(overlap),'p':float(hypergeom.sf(len(overlap)-1,len(background),len(genes),len(fg))),'overlap_genes':';'.join(sorted(overlap))})
out=pd.DataFrame(rows);out['global_functional_fdr']=multipletests(out.p,method='fdr_bh')[1];out=out.sort_values(['global_functional_fdr','p','term']);out.to_csv(R/'functional_GO_BP_all_tests.tsv.gz',sep='\t',index=False,compression='gzip')
top=out.groupby('foreground',sort=False).head(15);top.to_csv(R/'functional_GO_BP_top_terms.tsv',sep='\t',index=False)
for label in foregrounds:
 part=out[out.foreground.eq(label)];print(label,'genes',len(foregrounds[label]),'terms',len(part),'global FDR < .05',int(part.global_functional_fdr.lt(.05).sum()));print(part.head(8).to_string(index=False),flush=True)

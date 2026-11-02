from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import anndata as ad
from run_analysis import ROOT,RAW,P,R,annotation,gene_sets,design,ols_all

checks=[]
def check(name,passed):
 assert bool(passed),name
 checks.append({'check':name,'passed':True})

c=pd.read_csv(P/'GSE63514_metadata.tsv',sep='\t',index_col=0);primary=c.loc[c.primary_included]
check('Discovery primary contains 40 CIN3 and 20 cancer samples',primary.stage.value_counts().to_dict()=={'CIN3':40,'Cancer':20})
check('Discovery primary excludes all eight whole sections',primary.whole_section.eq(0).all() and c.whole_section.sum()==8)
c3=pd.read_csv(P/'GSE7803_metadata.tsv',sep='\t',index_col=0)
check('Three explicitly labeled technical/test arrays excluded',c3.excluded_test.sum()==3 and not (c3.excluded_test&c3.primary_included).any())
check('GSE7803 biological validation contains 7 HSIL and 21 SCC',c3.loc[c3.primary_included].stage.value_counts().to_dict()=={'Cancer':21,'HSIL':7})
d=pd.read_csv(R/'discovery_all_genes.tsv',sep='\t');check('Discovery BH corrects all 6685 primary tests',len(d)==6685 and np.allclose(multipletests(d.p,method='fdr_bh')[1],d.fdr))
excluded=set.union(*gene_sets().values());check('No tested gene belongs to any frozen proliferation marker set',not(set(d.gene)&excluded))
lock=json.loads((R/'candidate_lock.json').read_text());selected=d[d.fdr.lt(.05)&d.beta.abs().ge(.5)].head(5)
check('Frozen top-five candidates follow discovery-only ranking',[r['gene'] for r in lock['candidates']]==selected.gene.tolist())
check('Original plan hash retained',lock['plan_sha256']==hashlib.sha256((ROOT/'ANALYSIS_PLAN_v1.md').read_bytes()).hexdigest())
v=pd.read_csv(R/'external_validation.tsv',sep='\t')
for cohort,part in v.groupby('cohort'):
 check('Holm correction includes unmeasured targets in '+cohort,np.allclose(multipletests(part.p,method='holm')[1],part.holm_p))
check('Original screen is still reported as failed',not v.replicated.any())
ext=pd.read_csv(R/'exploratory_all_qualified_external.tsv',sep='\t')
for cohort,part in ext.groupby('cohort'):
 check('Exploratory FDR covers all 257 qualified targets in '+cohort,len(part)==257 and np.allclose(multipletests(part.p,method='fdr_bh')[1],part.extension_fdr))
check('Exploratory intersection exactly PLOD2 DEPTOR KLK11',set(pd.read_csv(R/'exploratory_three_cohort_intersection.tsv',sep='\t').gene)=={'PLOD2','DEPTOR','KLK11'})
rna=pd.read_csv(P/'GSE223804_metadata.tsv',sep='\t',index_col=0);rs=rna.loc[rna.primary_included]
check('RNA validation includes 4 HSIL and 9 invasive SCC',rs.stage.value_counts().to_dict()=={'Cancer':9,'HSIL':4})
check('RNA column CIN1 is CIN3 specimen and CIN5 is excluded low-grade CIN',rs.loc['CIN1','stage']=='HSIL' and not rna.loc['CIN5','primary_included'])
rl=json.loads((R/'RNAseq_candidate_lock.json').read_text());rr=pd.read_csv(R/'new_RNAseq_validation.tsv',sep='\t')
check('RNA tests exactly the three frozen extension genes',set(rr.gene)==set(rl['genes']))
check('New RNA Holm correction covers all three genes',np.allclose(multipletests(rr.p,method='holm')[1],rr.holm_p))
check('No RNA candidate is mislabeled as replicated',not rr.replicated.any())
check('Vectorized HC3 independently verified',json.loads((R/'numerical_verification.json').read_text())['passed'])
a=ad.read_h5ad(P/'GSE171894_selected_gene_localization.h5ad');check('Saved single-cell object reopens with 14452 cells from 4 donors',a.n_obs==14452 and a.obs.donor.nunique()==4)
check('All original and extension targets accounted for in localization',set(a.var_names)|set(json.loads((R/'singlecell_reference_provenance.json').read_text())['absent_genes'])==set(selected.gene)|set(rl['genes']))
records=[]
for f in sorted(RAW.glob('*.source.json')):
 r=json.loads(f.read_text());target=RAW/r['file'];check('Source checksum '+r['file'],hashlib.sha256(target.read_bytes()).hexdigest()==r['sha256']);records.append(r)
(R/'all_source_manifest.json').write_text(json.dumps(records,indent=2))
(R/'verification.json').write_text(json.dumps({'all_passed':True,'checks':checks},indent=2))
print('PASS',len(checks),'scientific, numerical, file-integrity and provenance checks')

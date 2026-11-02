from pathlib import Path
import hashlib,json,time
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT,R,dump

source=R/'genomewide_unadjusted_adjusted_coefficients.tsv.gz'
d=pd.read_csv(source,sep='\t');common=pd.read_csv(R/'genomewide_common_measured_universe.tsv',sep='\t').gene.tolist()
labels=['GSE63514','GSE138080','GSE7803'];ns=np.array([60,25,28]);weights=np.sqrt(ns)
p=d.pivot(index='gene',columns='cohort',values='p_adjusted').loc[common,labels]
b=d.pivot(index='gene',columns='cohort',values='beta_adjusted').loc[common,labels]
sb=d.pivot(index='gene',columns='cohort',values='standardized_beta_adjusted').loc[common,labels]
assert len(p)==6388 and p.notna().all().all()
z=np.sign(b)*norm.isf(np.clip(p.to_numpy(),1e-300,1)/2)
combined=np.asarray(z@weights)/np.sqrt(ns.sum());pp=2*norm.sf(np.abs(combined))
res=pd.DataFrame({'gene':common,'combined_z':combined,'combined_p':pp,'combined_fdr':multipletests(pp,method='fdr_bh')[1],'mean_absolute_standardized_beta':sb.abs().mean(axis=1).values,'same_direction_all_three':(np.sign(b).nunique(axis=1)==1).values,'direction':np.sign(b.iloc[:,0]).astype(int).values})
for label in labels:
 res['beta_'+label]=b[label].values;res['p_'+label]=p[label].values
res=res.sort_values(['combined_p','gene']);res.to_csv(R/'meta_discovery_all_genes.tsv',sep='\t',index=False)
qualified=res[res.same_direction_all_three&res.combined_fdr.lt(.05)&res.mean_absolute_standardized_beta.ge(.5)]
selected=qualified.head(10)
lockfile=R/'meta_discovery_candidate_lock.json'
record={'locked_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'new_cohort_target_models_not_yet_run':True,'format_preview_before_lock':'After the rule was written, the RNA count header and first three non-target rows were previewed during a failed lock-script run; no target comparisons or new-cohort outcome statistics were used. The bug fix only converts a pandas index to a numpy array.','analysis_role':'post hoc joint discovery in previously analyzed cohorts; new studies held out from ranking and model development','weights':dict(zip(labels,weights.tolist())),'genes_tested':len(res),'qualified':len(qualified),'candidates':selected.to_dict('records'),'plan_sha256':hashlib.sha256((ROOT/'META_DISCOVERY_EXTENSION_v1.md').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'validation':['GSE166466','GSE149763']}
if lockfile.exists():assert [x['gene'] for x in json.loads(lockfile.read_text())['candidates']]==selected.gene.tolist()
else:lockfile.write_text(json.dumps(record,indent=2))
print('Joint-discovery candidates fixed before new target-model testing; see timing note:')
print(selected.to_string(index=False));print('Total qualifying',len(qualified))

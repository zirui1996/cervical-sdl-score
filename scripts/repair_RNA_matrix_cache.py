"""Rebuild a truncated processed gzip from the unchanged source; verify effects."""
import hashlib,json
import numpy as np
import pandas as pd
from run_analysis import RAW,P,R,design,ols_all
from extension_utils import extra_scores
raw=pd.read_csv(RAW/'GSE223804_cpm_CIN_CR_Norm.txt.gz',sep='\t',index_col=0)
e=np.log2(raw+1).groupby(level=0).median();tmp=P/'GSE223804_log2_CPMplus1.atomic.tsv.gz';path=P/'GSE223804_log2_CPMplus1.tsv.gz'
oldhash=hashlib.sha256(path.read_bytes()).hexdigest()
e.to_csv(tmp,sep='\t',compression='gzip');check=pd.read_csv(tmp,sep='\t',index_col=0);assert check.shape==e.shape and np.allclose(e,check)
c=pd.read_csv(P/'GSE223804_metadata.tsv',sep='\t',index_col=0);cc=c.loc[c.primary_included];s,_=extra_scores(e)
old=pd.read_csv(R/'new_RNAseq_validation.tsv',sep='\t').set_index('gene');new=ols_all(e.loc[old.index],design(cc,s)).set_index('gene')
assert np.allclose(new.loc[old.index,'beta'],old.beta) and np.allclose(new.loc[old.index,'p'],old.p)
tmp.replace(path)
(R/'processed_RNA_cache_repair.json').write_text(json.dumps({'file':path.name,'reason':'prior processed gzip lacked an end-of-stream marker','source_unchanged':True,'original_three_gene_beta_and_p_reproduced':True,'reopened_shape':list(check.shape),'old_sha256':oldhash,'new_sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2))
print('Processed RNA matrix repaired atomically; all original three-gene effects and P values reproduced')

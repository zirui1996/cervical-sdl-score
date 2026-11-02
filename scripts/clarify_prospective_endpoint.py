import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from run_analysis import P,R
from extension_utils import frozen_targets,fit_targets,program_fit

c=pd.read_csv(P/'GSE75132_extension_metadata.tsv',sep='\t',index_col=0);e=pd.read_csv(P/'GSE75132_extension_gene_log_expression.tsv.gz',sep='\t',index_col=0);s=pd.read_csv(P/'GSE75132_extension_scores.tsv',sep='\t',index_col=0)
cross=pd.crosstab(c.category_code,c['disease state']);cross.to_csv(R/'GSE75132_category_disease_state_crosstab.tsv',sep='\t')
assert ((c.category_code==2)&c['disease state'].eq('moderate dysplasia')).sum()==4
old,meta,directions=frozen_targets();rows=[]
for label in ['GEO_progression_label_mixed_grade','Restricted_severe_CIS_cancer_annotations']:
 cc=c.copy()
 if label.startswith('Restricted'):cc['primary_included']=cc.primary_included&~cc['disease state'].eq('moderate dysplasia')
 tests=fit_targets('GSE75132',cc,e,s,list(directions),directions,True)+[program_fit('GSE75132',cc,s,True)]
 for record in tests:record['endpoint_definition']=label
 rows+=tests
result=pd.DataFrame(rows);result['holm_within_endpoint']=result.groupby('endpoint_definition').p.transform(lambda x:multipletests(x,method='holm')[1]);result['holm_all_28_endpoint_tests']=multipletests(result.p,method='holm')[1];result.to_csv(R/'prospective_endpoint_definition_sensitivity.tsv',sep='\t',index=False)
audit=json.loads((R/'GSE75132_label_audit.json').read_text());audit.update({'disease_state_in_progression_group':c.loc[c.category_code.eq(2),'disease state'].value_counts().to_dict(),'strict_CIN3plus_endpoint_not_established':True,'restricted_annotation_subset':'Exclude four moderate dysplasia annotations; 16 severe/CIS/cancer versus 10 none','limits':'No individual event times or diagnostic chronology. GEO category progression is mixed-grade; cannot equate every record with a confirmed strict CIN3+ event.'});(R/'GSE75132_label_audit.json').write_text(json.dumps(audit,indent=2))
primary=pd.read_csv(R/'prospective_endpoint_transfer.tsv',sep='\t');primary['endpoint_definition']='GEO_progression_label_mixed_grade';primary.to_csv(R/'prospective_endpoint_transfer.tsv',sep='\t',index=False)
print(cross.to_string());print(result[result.gene.eq('Joint_program')][['endpoint_definition','n','beta','p','ci_low','ci_high','holm_all_28_endpoint_tests']].to_string(index=False))

from pathlib import Path
import re,json
import pandas as pd
from run_analysis import ROOT,RAW,R,P

decisions={
'GSE63514':('joint discovery','40 CIN3 + 20 laser-captured cancers','Original primary specimens; eight whole sections remain sensitivity only'),
'GSE138080':('joint discovery','15 CIN2/3 + 10 SCC','Broader precursor category; older external data now used in joint discovery'),
'GSE7803':('joint discovery','7 HSIL + 21 SCC','Three test arrays excluded'),
'GSE223804':('prior RNA context','4 high-grade + 9 SCC','Initial frozen three-gene RNA validation; new program is post hoc context; ARCHS4 reprocessing is not a new cohort'),
'GSE166466':('new tissue validation','6 HPV16+ HSIL + 7 HPV16+ SCC; 7 normal context','Same Clariom D platform; no target/model tuning from these results'),
'GSE149763':('new RNA tissue validation','3 CIN3 + 3 cervical cancer; 3 cervicitis context','Very small; generic cervical-cancer label retained rather than assuming SCC from GEO alone'),
'GSE147009':('new RNA tissue validation','6 HPV16+ HSIL + 6 HPV16+ SCC; 12 non-lesional context','Uniformly processed host-gene counts obtained from ARCHS4, not author circRNA matrix'),
'GSE51993':('new precursor context','7 normal + 9 CIN1 + 8 CIN3','mRNA GPL10558 only; no invasive cancer'),
'GSE75132':('new prospective endpoint audit','10 HPV16-persistent without progression + 20 with progression label; 11 HPV-negative controls','Four progression records say moderate dysplasia; binary GEO labels are not verified strict CIN3+ events; no individual event times'),
'GSE208653':('new single-cell context','9 specimens: 2 HPV- normal + 2 HPV+ normal + 2 HSIL + 2 SCC + 1 ADC','Count-level QC; broad compartments and marker states; specimen-level descriptions, not cell-wise confirmation'),
'GSE171894':('existing single-cell reference','4 donors','Previously analyzed source context; not another precursor/invasion comparison'),
'GSE27678':('platform audit only','21 HSIL on GPL571; 28 SCC on GPL570','Primary HSIL/SCC comparison confounded with platform; normal bridge samples require extra assumptions, not primary confirmation'),
'GSE333892':('not analyzed','52 tissue specimens + organoid libraries in metadata','Public processed RDS about 21 GiB, already cell-cycle-regressed; 8 GiB runtime cannot support the intended unregressed-count analysis'),
'GSE120691':('excluded from tissue validation','11 cancers + 21 CIN + 19 controls','Peripheral blood cells; provided series files are differential-result tables; incompatible source'),
'GSE124288':('metadata only','2 preinvasive + 2 microinvasive','Insufficient size and potential source overlap with GSE223804 not resolved'),
'GSE64217':('metadata/small-data feasibility only','2 normal + 2 CIN + 2 SCC','Too small for powered validation'),
'GSE189452':('not included in main comparison','15 CIN3','HPV-integration comparison within CIN3; no invasive cancers'),
'GSE171597':('excluded from gene-expression analysis','miRNA LGCIN outcome study','miRNA rather than required mRNA assay'),
'GSE234041':('excluded from patient validation','18 cultured-cell libraries','Conditionally reprogrammed cells with estradiol perturbation; replicates are not independent primary-tissue cohorts'),
'GSE67522':('not included in invasion-boundary comparison','HPV-negative normal, HPV16+ nonmalignant and cancer','No explicitly comparable HSIL group'),
'GSE113942':('not included in invasion-boundary comparison','7 normal + 7 cancers','CIN validation was qPCR, not a deposited CIN RNA-seq matrix'),
'GSE145372':('not included in main SCC comparison','4 adenocarcinoma + 4 normal','Different histology and no HSIL'),
'GSE206224':('not included in main comparison','9 cancer/adjacent pairs','No HSIL; also circ/lncRNA/mRNA source'),
'GSE173112':('not included in main comparison','cancer/normal and cell-line perturbations','No HSIL; do not count knockdown replicates as patients'),
'GSE59322':('excluded from patient validation','sorted CIN612-9E cells','Cell-line sorting experiment'),
'GSE298989':('search false positive','THP-1/Candida experiment','Not cervical tissue lesions'),
'GSE274628':('search false positive','Fallopian-tube tissue','Different organ/question'),
'GSE264492':('excluded from patient validation','HeLa bacterial-infection time series','Cell-line perturbation'),
'GSE20499':('search false positive','decidual/endometrial NK cells','Different tissue/cell question'),
'GSE12852':('search false positive','pelvic-support tissue in prolapse','Different tissue/question'),
'GSE138081':('superseries duplicate','contains GSE138080','No additional independent patient cohort'),
'GSE27333':('superseries duplicate','contains GSE27678','No additional independent patient cohort'),
'GSE147011':('superseries duplicate','contains GSE147009','No additional independent patient cohort')}
meta={}
for f in RAW.glob('GSE*_metadata.txt'):
 text=f.read_text();match=re.search(r'!Series_geo_accession = (GSE\d+)',text)
 if match:
  accession='GSE'+str(int(match.group(1)[3:]));title=re.search(r'!Series_title = (.*)',text)
  meta[accession]=title.group(1).strip() if title else ''
rows=[]
for accession,(role,samples,reason) in decisions.items():rows.append({'accession':accession,'source_title':meta.get(accession,''),'role':role,'sample_context':samples,'decision_reason':reason,'source_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+accession})
pd.DataFrame(rows).to_csv(R/'cohort_eligibility_and_roles.tsv',sep='\t',index=False)
# Exact within-study profile duplication audit. Absence of exact duplicates is
# not proof that no patient appears in two separately deposited studies.
duplicate=[]
for accession in ['GSE166466','GSE149763','GSE147009','GSE51993','GSE75132']:
 e=pd.read_csv(P/f'{accession}_extension_gene_log_expression.tsv.gz',sep='\t',index_col=0)
 groups={}
 import hashlib
 for column in e:
  fingerprint=hashlib.sha256(e[column].to_numpy().tobytes()).hexdigest();groups.setdefault(fingerprint,[]).append(column)
 duplicate.append({'cohort':accession,'sample_columns':len(e.columns),'exact_duplicate_profile_sets':json.dumps([v for v in groups.values() if len(v)>1]),'duplicated_profile_columns':sum(len(v) for v in groups.values() if len(v)>1)})
pd.DataFrame(duplicate).to_csv(R/'additional_cohort_exact_profile_duplicates.tsv',sep='\t',index=False)
print('Cohort-role audit',len(rows),'records; new profile duplicates:',duplicate)

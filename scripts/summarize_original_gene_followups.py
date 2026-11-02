import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
from run_analysis import R
genes=['PLOD2','DEPTOR','KLK11']
new=pd.read_csv(R/'new_three_study_full_multiplicity.tsv',sep='\t');new=new[new.gene.isin(genes)][['cohort','gene','beta','p','direction_concordant']]
old=pd.read_csv(R/'new_RNAseq_validation.tsv',sep='\t');old=old[['gene','beta','p','direction_concordant']].copy();old['cohort']='GSE223804'
all=pd.concat([old,new],ignore_index=True);assert len(all)==12
all['holm_original_three_across_four_followups']=multipletests(all.p,method='holm')[1]
planned=all.cohort.isin(['GSE223804','GSE166466','GSE149763']);all['holm_planned_nine_test_subset']=np.nan;all.loc[planned,'holm_planned_nine_test_subset']=multipletests(all.loc[planned,'p'],method='holm')[1]
all.to_csv(R/'original_three_gene_all_followup_families.tsv',sep='\t',index=False)
print('All original three genes retained across four follow-up studies; no multiplicity-adjusted confirmation')

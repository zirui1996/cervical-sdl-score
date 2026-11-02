from pathlib import Path
import csv,gzip,json,re
import numpy as np
import pandas as pd
from run_analysis import ROOT,RAW,P,R,gene_sets,scores,design,ols_all,annotation,expression

def read_meta(filename):
 fields={};chars=[]
 with gzip.open(RAW/filename,'rt') as f:
  for line in f:
   if line.startswith('!series_matrix_table_begin'):break
   if not line.startswith('!Sample_'):continue
   row=next(csv.reader([line],delimiter='\t'));key=row[0][8:]
   if key=='characteristics_ch1':chars.append(row[1:])
   elif key not in fields:fields[key]=row[1:]
 c=pd.DataFrame(index=fields['geo_accession']);c.index.name='sample'
 for key in ['title','source_name_ch1','platform_id','description','data_processing']:
  if key in fields:c[key]=fields[key]
 for i,row in enumerate(chars):
  key=row[0].split(':',1)[0] if ':' in row[0] else 'characteristic_'+str(i)
  if key in c:key+='_'+str(i)
  c[key]=[x.split(':',1)[1].strip() if ':' in x else x for x in row]
 return c

def clariom_annotation():
 output=P/'GPL23126_unambiguous_probe_map.tsv'
 if output.exists():return pd.read_csv(output,sep='\t',index_col=0).gene
 rows=[];seen=False
 with open(RAW/'GPL23126_metadata.txt') as f:
  for line in f:
   if line.startswith('!platform_table_begin'):seen=True;header=next(f).rstrip('\r\n').split('\t');gi=header.index('gene_assignment');continue
   if line.startswith('!platform_table_end'):break
   if not seen:continue
   fields=line.rstrip('\r\n').split('\t')
   if len(fields)<=gi:continue
   symbols=set()
   for item in fields[gi].split(' /// '):
    vals=item.split(' // ')
    if len(vals)>=2:
     gene=vals[1].strip(' "')
     if re.fullmatch('[A-Za-z0-9_.-]+',gene) and gene!='---':symbols.add(gene)
   if len(symbols)==1:rows.append((fields[0],next(iter(symbols))))
 a=pd.DataFrame(rows,columns=['probe','gene']).drop_duplicates('probe');a.to_csv(output,sep='\t',index=False)
 return a.set_index('probe').gene

def common_universe():return pd.read_csv(R/'platform_common_gene_universe.tsv',sep='\t').gene.tolist()

def frozen_targets():
 old=json.loads((R/'RNAseq_candidate_lock.json').read_text())['genes']
 m=json.loads((R/'meta_discovery_candidate_lock.json').read_text())['candidates']
 directions={g:int(np.sign(pd.read_csv(R/'exploratory_three_cohort_intersection.tsv',sep='\t').set_index('gene').loc[g,'discovery_beta'])) for g in old}
 directions.update({x['gene']:int(x['direction']) for x in m})
 return old,[x['gene'] for x in m],directions

def extra_scores(e):
 s,cov=scores(e,common_universe(),gene_sets())
 rank=e.loc[e.index.intersection(common_universe())].rank(axis=0,pct=True)
 lock=json.loads((R/'joint_program_lock.json').read_text());valid=[x for x in lock['candidates'] if x['gene'] in rank.index]
 s['Joint_program']=np.nan
 if len(valid)>=lock['minimum_measured']:
  s['Joint_program']=sum(x['direction']*rank.loc[x['gene']] for x in valid)/len(valid)
 immune=['PTPRC','CD3D','CD3E','CD79A','MS4A1','LST1','TYROBP','FCER1G'];stroma=['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA']
 for name,genes in [('Immune_proxy',immune),('Stroma_proxy',stroma)]:
  have=rank.index.intersection(genes);s[name]=rank.loc[have].mean(axis=0)
  cov.loc[len(cov)]=[name,len(genes),len(have),';'.join(have),len(rank)]
 cov.loc[len(cov)]=['Joint_program',len(lock['candidates']),len(valid),';'.join(x['gene'] for x in valid),len(rank)]
 return s,cov

def progression_ratio_expression():
 # This source is a positive control-normalized RATIO matrix (median 1),
 # not an ordinary raw-intensity array. Preserve ratios below one.
 with gzip.open(RAW/'GSE75132_series_matrix.txt.gz','rt') as f:
  for line in f:
   if line.startswith('!series_matrix_table_begin'):break
  import io
  body=[]
  for line in f:
   if line.startswith('!series_matrix_table_end'):break
   body.append(line)
 e=pd.read_csv(io.StringIO(''.join(body)),sep='\t',index_col=0)
 assert (e.to_numpy()>0).all(),'Ratio zero/nonpositive values require source review'
 q=np.quantile(e.to_numpy(),[0,.01,.5,.99,1]);e=np.log2(e)
 ann=annotation('GPL570');symbols=ann.reindex(e.index.astype(str));e=e.loc[symbols.notna()].copy();e.index=symbols.dropna().to_numpy();e=e.groupby(level=0).median()
 record={'source_positive_ratio_quantiles':q.tolist(),'transformation':'log2 of every positive source ratio; no floor at 1, no pseudocount','rationale':'Author data were divided by HPV-negative reference medians; values below one carry down-regulation information.','genes':len(e),'samples':e.shape[1]}
 (R/'GSE75132_ratio_scale_audit.json').write_text(json.dumps(record,indent=2));(R/'GSE75132_processing.json').write_text(json.dumps(record,indent=2))
 tmp=P/'GSE75132_gene_log_expression.atomic.tsv.gz';e.to_csv(tmp,sep='\t',compression='gzip');assert pd.read_csv(tmp,sep='\t',index_col=0).shape==e.shape;tmp.replace(P/'GSE75132_gene_log_expression.tsv.gz');return e

def prepare_new_bulk():
 cohorts={}
 # The missing standard annot.gz is replaced by a streamed parse of GEO's
 # complete, downloaded native platform table, preserving unambiguous symbols.
 ann=clariom_annotation();e=expression('GSE166466',ann);c=read_meta('GSE166466_series_matrix.txt.gz')
 c['stage']=np.select([c.title.str.contains('HSIL'),c.title.str.contains('squamous cell carcinoma')],['HSIL','Cancer'],default='Normal');c['group']=c.stage.eq('Cancer').astype(int);c['primary_included']=c.stage.isin(['HSIL','Cancer']);cohorts['GSE166466']=(c,e)
 reads=pd.read_csv(RAW/'GSE149763_expressed_gene_reads.txt.gz',sep='\t');rpkm=pd.read_csv(RAW/'GSE149763_expressed_gene_RPKM.txt.gz',sep='\t');c=read_meta('GSE149763_series_matrix.txt.gz')
 samplecols=c.title.tolist();ct=reads[samplecols].apply(pd.to_numeric)
 assert (ct.to_numpy()>=0).all() and np.allclose(ct,ct.round()),'Read matrix units need review'
 library=ct.sum(axis=0);symbols=reads.Symbol.astype(str);ok=symbols.str.fullmatch('[A-Za-z0-9_.-]+')&~symbols.eq('---')
 ct=ct.loc[ok];ct.index=symbols[ok];ct=ct.groupby(level=0).sum();ct.columns=c.index
 e=np.log2(ct.div(library.to_numpy(),axis=1)*1e6+1);ct.to_csv(P/'GSE149763_mapped_gene_reads.tsv.gz',sep='\t',compression='gzip')
 r=rpkm[samplecols].apply(pd.to_numeric);r.index=rpkm.Symbol.astype(str);r=r.loc[r.index.str.fullmatch('[A-Za-z0-9_.-]+')&(r.index!='---')].groupby(level=0).median();r.columns=c.index
 np.log2(r+1).to_csv(P/'GSE149763_author_log_RPKM.tsv.gz',sep='\t',compression='gzip')
 (R/'GSE149763_count_scale_audit.json').write_text(json.dumps({'nonnegative_integer_source_reads':True,'library_sums_all_source_genes':dict(zip(c.index,library.astype(int).tolist())),'mapping':'sum counts for a unique source symbol; no ambiguous symbol mappings','normalization':'log2(1+CPM) with full source gene-read column totals; not assumed to reproduce author edgeR pipeline'},indent=2))
 c['stage']=np.select([c.title.str.startswith('CINIII'),c.title.str.startswith('IIa')],['HSIL','Cancer'],default='Cervicitis');c['group']=c.stage.eq('Cancer').astype(int);c['primary_included']=c.stage.isin(['HSIL','Cancer']);cohorts['GSE149763']=(c,e)
 # Existing expression() accepts the platform-suffixed source basename.
 en=expression('GSE51993-GPL10558',annotation('GPL10558'));c=read_meta('GSE51993-GPL10558_series_matrix.txt.gz')
 st=c['disease stage'];c['stage']=np.select([st.str.contains('III'),st.str.contains('CIN')],['HSIL','CIN1'],default='Normal');c['group']=c.stage.eq('HSIL').astype(int);c['primary_included']=c.stage.isin(['Normal','HSIL']);cohorts['GSE51993']=(c,en)
 e=progression_ratio_expression();c=read_meta('GSE75132_series_matrix.txt.gz');key=[x for x in c.columns if x.startswith('category (')][0]
 c['category_code']=c[key].astype(int);c['stage']=c.category_code.map({0:'HPV_negative',1:'Persistent_no_progression',2:'Persistent_progression'});c['group']=c.category_code.eq(2).astype(int);c['primary_included']=c.category_code.gt(0)
 titlecode=np.select([c.title.str.contains('without progression'),c.title.str.contains('_with progression')],[1,2],default=0)
 assert np.array_equal(titlecode,c.category_code),'GEO title/category disagreement'
 (R/'GSE75132_label_audit.json').write_text(json.dumps({'GEO_category_counts':c.category_code.value_counts().sort_index().to_dict(),'title_category_agreement':True,'progression_field_definition':key,'primary_source_paper':'PMID 27508094; normal cytology at baseline; paper follow-up endpoint CIN3+','secondary_count_disagreement':'Some secondary sources give different 20/10 counts. Preserve GEO source labels; their endpoint is not automatically strict CIN3+.','disease_state_in_progression_group':c.loc[c.category_code.eq(2),'disease state'].value_counts().to_dict(),'strict_CIN3plus_endpoint_not_established':True,'limits':'No individual times or diagnostic chronology; four category=2 records have moderate dysplasia annotation. Do not equate every progression label with a confirmed strict CIN3+ event.'},indent=2));cohorts['GSE75132']=(c,e)
 for label,(c,e) in cohorts.items():
  s,cov=extra_scores(e);c.to_csv(P/f'{label}_extension_metadata.tsv',sep='\t');e.to_csv(P/f'{label}_extension_gene_log_expression.tsv.gz',sep='\t',compression='gzip');s.to_csv(P/f'{label}_extension_scores.tsv',sep='\t');cov.to_csv(R/f'{label}_extension_score_coverage.tsv',sep='\t',index=False)
  cohorts[label]=(c,e,s)
 return cohorts

def fit_targets(label,c,e,s,targets,directions,adjusted=True):
 cc=c.loc[c.primary_included].copy();X=design(cc,s) if adjusted else pd.DataFrame({'const':1.,'group':cc.group.astype(float)},index=cc.index)
 rows=[]
 for gene in targets:
  record={'cohort':label,'gene':gene,'model':'proliferation_adjusted' if adjusted else 'unadjusted','n':len(cc),'precursor_or_comparator':int(cc.group.eq(0).sum()),'invasive_or_event':int(cc.group.sum()),'expected_direction':directions[gene]}
  if gene not in e.index or e.loc[gene,cc.index].isna().any() or e.loc[gene,cc.index].std()==0:
   record.update({'p':1.,'status':'unmeasured_or_nonestimable','direction_concordant':False})
  else:
   rr=ols_all(e.loc[[gene]],X).iloc[0].to_dict();rr.pop('fdr');record.update(rr);record.update({'status':'tested','direction_concordant':bool(np.sign(rr['beta'])==directions[gene])})
  rows.append(record)
 return rows

def program_fit(label,c,s,adjusted=True):
 cc=c.loc[c.primary_included].copy();v=s.loc[cc.index,'Joint_program'];record={'cohort':label,'gene':'Joint_program','model':'proliferation_adjusted' if adjusted else 'unadjusted','n':len(cc),'expected_direction':1}
 if v.isna().any() or v.std()==0:return {**record,'p':1.,'status':'unscorable','direction_concordant':False}
 ev=pd.DataFrame([(v-v.mean())/v.std(ddof=1)],index=['Joint_program'])
 XX=design(cc,s) if adjusted else pd.DataFrame({'const':1.,'group':cc.group.astype(float)},index=cc.index)
 rr=ols_all(ev,XX).iloc[0].to_dict();rr.pop('fdr');return {**record,**rr,'status':'tested','direction_concordant':bool(rr['beta']>0)}

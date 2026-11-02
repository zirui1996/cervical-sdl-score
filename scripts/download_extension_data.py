"""Download additional cohort matrices and exact public single-cell sample files."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,sys
from download_data import download,RAW,ROOT

SOURCES={
 'GSE166466_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166466/matrix/GSE166466_series_matrix.txt.gz',
 # This platform has no standard annot.gz. Its native table contains the
 # transcript-to-symbol assignments and is streamed by the analysis code.
 'GPL23126_metadata.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23126&targ=self&form=text&view=full',
 'GPL10558.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz',
 'GSE51993-GPL10558_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51993/matrix/GSE51993-GPL10558_series_matrix.txt.gz',
 'GSE75132_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75132/matrix/GSE75132_series_matrix.txt.gz',
 'GSE149763_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149763/matrix/GSE149763_series_matrix.txt.gz',
 'GSE149763_expressed_gene_reads.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149763/suppl/GSE149763_expressed_gene_reads.txt.gz',
 'GSE149763_expressed_gene_RPKM.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149763/suppl/GSE149763_expressed_gene_RPKM.txt.gz',
 'GSE27678-GPL570_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27678/matrix/GSE27678-GPL570_series_matrix.txt.gz',
 'GSE27678-GPL571_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27678/matrix/GSE27678-GPL571_series_matrix.txt.gz'}
SOURCES['GSE208653_filelist.txt']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208653/suppl/filelist.txt'
SOURCES['GSE147009_metadata.txt']='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147009&targ=self&form=text&view=full'

def sc_sources():
 # Names are copied from the downloaded official GEO filelist, not inferred
 # histology. Preserve ADC_6 as adenocarcinoma, excluded from SCC contrasts.
 out={}
 for line in (RAW/'GSE208653_filelist.txt').read_text().splitlines():
  cols=line.split('\t')
  if len(cols)>1 and cols[0]=='File':
   name=cols[1];sample=name.split('_')[0]
   out[name]=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6360nnn/{sample}/suppl/{name}'
 return out

def main():
 src=sc_sources() if '--singlecell' in sys.argv else SOURCES
 results=[]
 with ThreadPoolExecutor(max_workers=4) as ex:
  pending={ex.submit(download,item):item for item in src.items()}
  for future in as_completed(pending):
   name,url=pending[future]
   try:results.append({'status':'downloaded',**future.result()})
   except Exception as e:
    record={'status':'failed','file':name,'url':url,'error':repr(e)};results.append(record);print(record,flush=True)
 target='extension_singlecell_download_manifest.json' if '--singlecell' in sys.argv else 'extension_bulk_download_manifest.json'
 (ROOT/'results'/target).write_text(json.dumps(results,indent=2))
 if any(r['status']=='failed' for r in results):raise SystemExit('Some sources failed; inspect manifest without substituting data.')

if __name__=='__main__':main()

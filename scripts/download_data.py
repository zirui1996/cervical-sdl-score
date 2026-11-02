from pathlib import Path
import concurrent.futures,hashlib,json,time,urllib.request
ROOT=Path(__file__).resolve().parents[1];RAW=ROOT/'data/raw';RAW.mkdir(parents=True,exist_ok=True)
SOURCES={
 'GSE63514_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63514/matrix/GSE63514_series_matrix.txt.gz',
 'GSE138080_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138080/matrix/GSE138080_series_matrix.txt.gz',
 'GSE7803_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7803/matrix/GSE7803_series_matrix.txt.gz',
 'GPL570.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz',
 'GPL4133.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4133/annot/GPL4133.annot.gz',
 'GPL96.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz',
 'h.all.v2025.1.Hs.symbols.gmt':'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt',
 'GSE223804_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223804/matrix/GSE223804_series_matrix.txt.gz',
 'GSE223804_cpm_CIN_CR_Norm.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223804/suppl/GSE223804_cpm_CIN_CR_Norm.txt.gz',
 'GSE64217_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64217/matrix/GSE64217_series_matrix.txt.gz'}
def download(item):
 name,url=item;f=RAW/name
 if not f.exists():
  with urllib.request.urlopen(url,timeout=120) as response,open(str(f)+'.part','wb') as out:
   while b:=response.read(1024*1024):out.write(b)
  Path(str(f)+'.part').rename(f)
 checksum=hashlib.sha256(f.read_bytes()).hexdigest()
 record={'file':name,'url':url,'bytes':f.stat().st_size,'sha256':checksum,'verified_or_retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
 (RAW/(name+'.source.json')).write_text(json.dumps(record,indent=2));print(name,record['bytes'],flush=True);return record
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(download,SOURCES.items()))
 (ROOT/'results/download_manifest.json').write_text(json.dumps(records,indent=2))

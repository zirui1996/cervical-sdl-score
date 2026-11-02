from pathlib import Path
import hashlib,json,time,urllib.request,sys
from download_data import RAW,ROOT

BASE='https://maayanlab.cloud/sigpy'
def save(name,data,url,method='GET',payload=None):
 f=RAW/name;f.write_bytes(data);r={'file':name,'url':url,'method':method,'payload':payload,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())};Path(str(f)+'.source.json').write_text(json.dumps(r,indent=2));return data

def main():
 accession=next((x for x in sys.argv[1:] if x.startswith('GSE')),'GSE147009')
 dest=RAW/(accession+'_ARCHS4_expression.zip')
 if dest.exists():
  import zipfile
  with zipfile.ZipFile(dest) as z:assert z.testzip() is None
  print('Existing archived expression reused; no network task needed',accession,dest.stat().st_size);return
 meta_file=RAW/(accession+'_ARCHS4_metadata_search.json')
 if not meta_file.exists():
  url=BASE+'/meta/quicksearch?query='+accession+'&species=human'
  with urllib.request.urlopen(url,timeout=90) as r:save(meta_file.name,r.read(),url)
 samples=json.loads(meta_file.read_text())['samples'];assert samples
 expected={'GSE147009':24,'GSE149763':9,'GSE223804':16};assert len(samples)==expected[accession],(accession,len(samples))
 taskfile=RAW/(accession+'_ARCHS4_expression_task.json');payload={'gsm_ids':samples,'species':'human'}
 if not taskfile.exists():
  url=BASE+'/data/samples';req=urllib.request.Request(url,data=json.dumps(payload).encode(),headers={'Content-Type':'application/json'},method='POST')
  with urllib.request.urlopen(req,timeout=120) as r:data=r.read()
  save(taskfile.name,data,url,'POST',payload)
 result=json.loads(taskfile.read_text());print('Task response',result,flush=True)
 task=result.get('task_id',result.get('id'));assert task,'Official API response has no task id'
 url=BASE+'/data/samples/status/'+task
 with urllib.request.urlopen(url,timeout=60) as r:status=json.loads(save(accession+'_ARCHS4_expression_status.json',r.read(),url))
 print('Task status',status,flush=True)
 if str(status.get('status','')).lower() not in ['complete','completed','success','finished','done']:
  print('Pending: run this script again to resume the same task; do not create a duplicate request.');return
 if not dest.exists():
  url=BASE+'/data/samples/download/'+task
  with urllib.request.urlopen(url,timeout=180) as r:data=r.read()
  save(dest.name,data,url)
 print('Saved expression archive',dest.stat().st_size)

if __name__=='__main__':main()

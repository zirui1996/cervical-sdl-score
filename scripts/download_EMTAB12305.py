"""Fetch official sparse count files, retaining byte sizes and SHA-256 provenance."""
import concurrent.futures,gzip,hashlib,json,time,urllib.request
from pathlib import Path
from run_analysis import RAW,R

DEST=RAW/'E-MTAB-12305';DEST.mkdir(exist_ok=True)
BASE='https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/305/E-MTAB-12305/Files/'

def download(name,expected=None):
    p=DEST/name
    if p.exists() and (expected is None or p.stat().st_size==expected):return {'file':name,'status':'cached','bytes':p.stat().st_size}
    errors=[]
    for attempt in range(3):
        for root in [BASE,'https://www.ebi.ac.uk/biostudies/files/E-MTAB-12305/']:
            url=root+name;tmp=p.with_name(p.name+'.partial')
            try:
                req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
                digest=hashlib.sha256();size=0
                with urllib.request.urlopen(req,timeout=30) as response,tmp.open('wb') as out:
                    while chunk:=response.read(1024*1024):out.write(chunk);digest.update(chunk);size+=len(chunk)
                if expected is not None:assert size==expected,(name,size,expected)
                if name.endswith('.gz'):
                    with gzip.open(tmp,'rb') as f:
                        while f.read(4*1024*1024):pass
                tmp.replace(p)
                record={'file':name,'url':url,'bytes':size,'sha256':digest.hexdigest(),'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'official_EBI_source':True}
                Path(str(p)+'.source.json').write_text(json.dumps(record,indent=2))
                print('Downloaded',name,size,flush=True);return {'file':name,'status':'downloaded','bytes':size}
            except Exception as e:
                errors.append({'url':url,'error':str(e)})
                if tmp.exists():tmp.unlink()
    print('Unavailable',name,errors[-1],flush=True)
    return {'file':name,'status':'failed','errors':errors}

def main():
    meta=json.loads((DEST/'metadata.json').read_text());files={}
    def walk(x):
        if isinstance(x,dict):
            if x.get('type')=='file' and 'path' in x:files[x['path']]=x['size']
            for v in x.values():walk(v)
        elif isinstance(x,list):
            for v in x:walk(v)
    walk(meta)
    records=[]
    for name in ['E-MTAB-12305.sdrf.txt','E-MTAB-12305.idf.txt']:
        records.append(download(name,files[name]))
    # Prioritize HSIL+tumor files, then adjacent context. Never fetch redundant
    # dense CSVs or lymph-node cells for this cervical-tissue question.
    jobs=[(sample+suffix,files[sample+suffix]) for sample in ['H1','H2','T1','T2','T3','T4','N1','N2','N3'] for suffix in ['features.tsv.gz','barcodes.tsv.gz','matrix.mtx.gz']]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for record in pool.map(lambda x:download(*x),jobs):records.append(record)
    (R/'EMTAB12305_download_audit.json').write_text(json.dumps(records,indent=2))
    print('Download audit',len(records),'files; failed',sum(x['status']=='failed' for x in records),flush=True)

if __name__=='__main__':main()

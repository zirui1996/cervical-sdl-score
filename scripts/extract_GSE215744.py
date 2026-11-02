import tarfile,openpyxl,io,csv,gzip
from pathlib import Path
root=Path(__file__).resolve().parents[1]
with tarfile.open(root/'data/raw/GSE215744_RAW.tar') as t:
 for m in t.getmembers():
  w=openpyxl.load_workbook(io.BytesIO(t.extractfile(m).read()),read_only=True,data_only=True)
  rows=w.active.iter_rows(values_only=True);header=next(rows)
  gi=header.index('Symbol');ci=next(i for i,h in enumerate(header) if str(h).endswith(':ReadCount'));ti=next(i for i,h in enumerate(header) if str(h).endswith(':RPKM'))
  with gzip.open(root/'data/processed'/m.name.replace('.xlsx','.tsv.gz'),'wt') as f:
   writer=csv.writer(f,delimiter='\t');writer.writerow(['gene','count','RPKM'])
   for row in rows:writer.writerow([row[gi],row[ci],row[ti]])
  w.close()

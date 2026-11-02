"""Retain or write the prespecified three-cohort gene lock before new RNA expression testing."""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';file=R/'RNAseq_candidate_lock.json'
genes=json.loads((R/'exploratory_extension_summary.json').read_text())['intersection_genes'][:10]
if file.exists():
 old=json.loads(file.read_text());assert old['genes']==genes,'New-cohort lock mismatch';print('Existing RNA candidate lock preserved.')
else:
 record={'locked_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'genes':genes,'selection':'all three-cohort intersection genes; no literature filter','validation':'GSE223804','plan_sha256':hashlib.sha256((ROOT/'NEW_COHORT_VALIDATION_RULES.md').read_bytes()).hexdigest(),'expression_not_yet_opened':True}
 file.write_text(json.dumps(record,indent=2));print('New RNA candidate lock written.')

"""Annotation only: official old/new platform bridge by unchanged oligo ID."""
import io
import pandas as pd
from run_analysis import RAW,P

def read_table(platform):
 text=(RAW/f'{platform}_metadata.txt').read_text().split('!platform_table_begin\n')[1].split('!platform_table_end')[0]
 return pd.read_csv(io.StringIO(text),sep='\t',dtype=str)

def main():
 a=read_table('GPL6699');b=read_table('GPL27630')
 b=b[b.Gene_Symbol.str.fullmatch('[A-Za-z0-9_.-]+',na=False)&b.Gene_Symbol.ne('---')]
 grouped=b.groupby('SPOT_ID').Gene_Symbol.agg(lambda x:sorted(set(x)))
 grouped=grouped[grouped.map(len).eq(1)].map(lambda x:x[0])
 m=a[['ID','Reporter Name']].assign(gene=a['Reporter Name'].map(grouped)).dropna()
 m[['ID','gene']].to_csv(P/'GPL6699_probe_map_via_official_alternative.tsv',sep='\t',index=False)

if __name__=='__main__':main()

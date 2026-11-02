import json
import numpy as np
import pandas as pd
from run_analysis import R,P,design,ols_all
from extension_utils import extra_scores,common_universe
from analyze_stage_and_composition import load_all

cohorts=load_all();rows=[];diag=[];lock=json.loads((R/'joint_program_lock.json').read_text());targets=[x['gene'] for x in lock['candidates']];directions={x['gene']:x['direction'] for x in lock['candidates']}
for label in ['GSE63514','GSE138080','GSE7803','GSE166466','GSE147009','GSE149763','GSE223804']:
 c,e,s=cohorts[label];cc=c.loc[c.primary_included].copy();v=s.loc[cc.index,'Joint_program'];Y=pd.DataFrame([(v-v.mean())/v.std(ddof=1)],index=['Joint_program'])
 def add(name,csub,X,yy=None):
  rr=ols_all(Y if yy is None else yy,X).iloc[0].to_dict();rr.update({'cohort':label,'sensitivity':name,'n':len(csub),'direction_concordant':bool(rr['beta']>0)});rows.append(rr)
 for kind in ['Composite','E2F','G2M','Fixed']:add(kind,cc,design(cc,s,kind))
 add('quadratic',cc,design(cc,s,quadratic=True))
 lo=max(s.loc[cc.index[cc.group.eq(k)],'Composite'].min() for k in [0,1]);hi=min(s.loc[cc.index[cc.group.eq(k)],'Composite'].max() for k in [0,1]);sub=cc.loc[s.loc[cc.index,'Composite'].between(lo,hi)]
 if sub.group.nunique()==2 and sub.group.value_counts().min()>=5:add('common_support',sub,design(sub,s))
 for sample in cc.index:
  sub=cc.drop(index=sample);add('leave_sample_'+sample,sub,design(sub,s))
 rank=e.loc[e.index.intersection(common_universe()),cc.index].rank(axis=0,pct=True)
 for target in targets:
  present=[g for g in targets if g!=target and g in rank.index]
  vv=sum(directions[g]*rank.loc[g] for g in present)/len(present);YY=pd.DataFrame([(vv-vv.mean())/vv.std(ddof=1)],index=['Joint_program']);add('leave_target_'+target,cc,design(cc,s),YY)
 # These sample-axis aggregation checks are post hoc and never replace the
 # frozen primary within-sample-rank program.
 vals=e.loc[[g for g in targets if g in e.index],cc.index].copy();sign=pd.Series(directions).reindex(vals.index)
 z=vals.sub(vals.mean(axis=1),axis=0).div(vals.std(axis=1).replace(0,np.nan),axis=0)
 for name,mat in [('post_hoc_gene_z_mean',z),('post_hoc_gene_sample_percentile_mean',vals.rank(axis=1,pct=True))]:
  vv=mat.mul(sign,axis=0).mean(axis=0);YY=pd.DataFrame([(vv-vv.mean())/vv.std(ddof=1)],index=['Joint_program']);add(name,cc,design(cc,s),YY)
 X=design(cc,s);rho=np.corrcoef(cc.group,s.loc[cc.index,'Composite'])[0,1]
 diag.append({'cohort':label,'n':len(cc),'group_proliferation_r':rho,'group_proliferation_VIF':1/(1-rho**2),'condition_number':np.linalg.cond(X),'precursor_common_support':int(((cc.group==0)&s.loc[cc.index,'Composite'].between(lo,hi)).sum()),'cancer_common_support':int(((cc.group==1)&s.loc[cc.index,'Composite'].between(lo,hi)).sum())})
d=pd.DataFrame(rows);d.to_csv(R/'joint_program_supportive_robustness.tsv',sep='\t',index=False);pd.DataFrame(diag).to_csv(R/'joint_program_design_diagnostics.tsv',sep='\t',index=False)
summary=[]
for label,part in d.groupby('cohort'):
 for name,mask in [('leave_one_sample',part.sensitivity.str.startswith('leave_sample_')),('leave_one_target',part.sensitivity.str.startswith('leave_target_')),('other',~part.sensitivity.str.startswith(('leave_sample_','leave_target_')))]:
  v=part.loc[mask];summary.append({'cohort':label,'check_family':name,'fits':len(v),'direction_positive_fraction':v.direction_concordant.mean(),'nominal_p_lt_05_fraction':v.p.lt(.05).mean(),'beta_min':v.beta.min(),'beta_max':v.beta.max()})
summary=pd.DataFrame(summary);summary.to_csv(R/'joint_program_robustness_summary.tsv',sep='\t',index=False);print(summary.to_string(index=False));print('\nNew validation supportive models:\n',d[d.cohort.eq('GSE166466')&~d.sensitivity.str.startswith('leave_sample_')][['sensitivity','beta','p','direction_concordant']].to_string(index=False))

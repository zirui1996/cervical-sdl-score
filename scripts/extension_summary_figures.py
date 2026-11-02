import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import R,F

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
primary=pd.read_csv(R/'joint_program_all_tissue_cohorts.tsv',sep='\t').set_index('cohort').loc['GSE166466']
comp=pd.read_csv(R/'composition_proxy_sensitivity.tsv',sep='\t');comp=comp[comp.cohort.eq('GSE166466')&comp.gene.eq('Joint_program')].iloc[0]
epi=pd.read_csv(R/'epithelial_context_proxy_sensitivity.tsv',sep='\t');epi=epi[epi.cohort.eq('GSE166466')&epi.gene.eq('Joint_program')]
records=[('Composite only',primary.beta,primary.ci_low,primary.ci_high),('Plus immune and stroma proxies',comp.beta,comp.ci_low,comp.ci_high)]
for _,row in epi.iterrows():records.append(('Plus '+row.proxy.replace('_',' '),row.beta,row.ci_low,row.ci_high))
fig,ax=plt.subplots(figsize=(8.5,4.4))
for i,(name,b,lo,hi) in enumerate(records):ax.errorbar(b,i,xerr=[[b-lo],[hi-b]],fmt='o',capsize=3,color='#3c657c')
ax.axvline(0,color='grey',ls='--');ax.set(yticks=range(len(records)),yticklabels=[x[0] for x in records],xlabel='Adjusted standardized program coefficient (95% HC3 CI)',title='GSE166466: post hoc biological-context sensitivities');ax.invert_yaxis();fig.tight_layout();fig.savefig(F/'14_program_epithelial_context.png',dpi=300);fig.savefig(F/'14_program_epithelial_context.pdf');plt.close(fig)
d=pd.read_csv(R/'GSE208653_epithelial_marker_states.tsv',sep='\t');stateorder=['Glandular_secretory','Squamous_basal','Squamous_differentiated','Mixed_unspecified'];names=['Secretory markers','Basal squamous markers','Differentiated squamous markers','Mixed or unspecified'];colors=['#7ba2b9','#ad6c50','#c5aa6b','#b9bec2']
fig,axs=plt.subplots(1,2,figsize=(12.2,4.9));p=d.pivot(index='donor',columns='marker_state',values='fraction_of_epithelial').fillna(0).reindex(columns=stateorder);bottom=np.zeros(len(p));labels=d.drop_duplicates('donor').set_index('donor').reindex(p.index)
for state,name,col in zip(stateorder,names,colors):axs[0].bar(np.arange(len(p)),p[state],bottom=bottom,color=col,label=name);bottom+=p[state].values
axs[0].set(xticks=range(len(p)),xticklabels=['N-1','N-2','N+1','N+2','HSIL1','HSIL2','SCC1','SCC2','ADC'],ylabel='Fraction of marker-defined epithelial cells',ylim=(0,1),title='A  Epithelial marker states vary by specimen');axs[0].tick_params(axis='x',labelsize=8,rotation=25)
sub=d[d.marker_state.eq('Squamous_basal')&d.cells.ge(30)];groups=['HPV_negative_normal','HPV_positive_normal','HSIL','SCC','ADC']
for i,group in enumerate(groups):
 v=sub[sub.stage.eq(group)];xx=np.linspace(i-.09,i+.09,len(v)) if len(v)>1 else [i]
 if len(v):
  axs[1].scatter(xx,v.Joint_program,s=40,color='#3d647d')
  for x,(_,row) in zip(xx,v.iterrows()):axs[1].annotate(str(row.cells),(x,row.Joint_program),xytext=(0,5),textcoords='offset points',fontsize=7,ha='center')
axs[1].set(xticks=range(5),xticklabels=['HPV- normal','HPV+ normal','HSIL','SCC','ADC'],ylabel='Frozen program signed mean fractional rank',title='B  Basal-marker state: one point per specimen');axs[1].tick_params(axis='x',rotation=25,labelsize=8)
handles,legend_labels=axs[0].get_legend_handles_labels();fig.legend(handles,legend_labels,loc='lower center',bbox_to_anchor=(.5,.01),ncol=2,frameon=False,fontsize=8)
fig.tight_layout(rect=(0,.14,1,1));fig.savefig(F/'15_singlecell_epithelial_states.png',dpi=300);fig.savefig(F/'15_singlecell_epithelial_states.pdf');plt.close(fig)
print('Generated epithelial-context and specimen-state figures')

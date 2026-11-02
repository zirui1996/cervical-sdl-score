import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import R,F

def main():
    out=pd.read_csv(R/'EMTAB12305_specimen_compartment_programs.tsv',sep='\t')
    fig,axs=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
    for ax,comp in zip(axs,['All cells','Epithelial']):
        for j,(group,color) in enumerate({'Adjacent_normal':'#74939b','HSIL':'#c8a04e','SCC':'#a95e51'}.items()):
            part=out[out.compartment.eq(comp)&out.stage.eq(group)&out.cells.ge(30)]
            x=np.linspace(j-.13,j+.13,len(part));ax.scatter(x,part.Joint_program,color=color,s=42)
            for xx,(_,row) in zip(x,part.iterrows()):ax.annotate(row.specimen,(xx,row.Joint_program),xytext=(3,4),textcoords='offset points',fontsize=7)
        ax.set(xticks=range(3),xticklabels=['Adjacent normal','HSIL','SCC'],ylabel='Frozen signed mean fractional rank',title=comp+' | one point per specimen',xlim=(-.35,2.45))
    fig.savefig(F/'18_EMTAB12305_fixed_program.png',dpi=300);fig.savefig(F/'18_EMTAB12305_fixed_program.pdf');plt.close(fig)
    old=pd.read_csv(R/'GSE208653_epithelial_marker_states.tsv',sep='\t').rename(columns={'donor':'specimen'});old=old[old.stage.isin(['HSIL','SCC'])].copy();old['source']='GSE208653';old['label']=old.specimen.str[-3:]
    new=pd.read_csv(R/'EMTAB12305_epithelial_marker_states.tsv',sep='\t');new['source']='E-MTAB-12305';new['label']=new.specimen
    combined=pd.concat([old,new],ignore_index=True);combined.to_csv(R/'two_source_epithelial_marker_state_comparison.tsv',sep='\t',index=False)
    names=['GSM6360684','GSM6360685','GSM6360686','GSM6360687','H1','H2','T1','T2','T3','T4']
    states=['Glandular_secretory','Squamous_basal','Squamous_differentiated','Mixed_unspecified'];colors=['#7d9eaa','#c59b58','#a66d60','#c4c9cc']
    fig,axs=plt.subplots(1,2,figsize=(12.2,4.8),layout='constrained',gridspec_kw={'width_ratios':[1.35,1]})
    bottom=np.zeros(len(names))
    for state,color in zip(states,colors):
        vals=combined[combined.marker_state.eq(state)].set_index('specimen').fraction_of_epithelial.reindex(names).fillna(0).to_numpy()
        axs[0].bar(range(len(names)),vals,bottom=bottom,color=color,label=state.replace('_',' '));bottom+=vals
    axs[0].axvline(3.5,color='grey',ls='--',lw=.8)
    axs[0].set(xticks=range(len(names)),xticklabels=['684 H','685 H','686 S','687 S','H1','H2','T1','T2','T3','T4'],ylabel='Fraction of marker-assigned epithelial cells',title='A  Fixed marker states across two sources')
    axs[0].tick_params(axis='x',rotation=50);axs[0].legend(frameon=False,fontsize=7,loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2)
    basal=combined[combined.marker_state.eq('Squamous_basal')&combined.cells.ge(30)]
    for j,(source,group) in enumerate([('GSE208653','HSIL'),('GSE208653','SCC'),('E-MTAB-12305','HSIL'),('E-MTAB-12305','SCC')]):
        d=basal[(basal.source==source)&(basal.stage==group)];x=np.linspace(j-.12,j+.12,len(d))
        axs[1].scatter(x,d.Joint_program,color='#c59b58' if group=='HSIL' else '#a66d60',s=38)
        for xx,(_,row) in zip(x,d.iterrows()):axs[1].annotate(row.label,(xx,row.Joint_program),xytext=(3,3),textcoords='offset points',fontsize=7)
    axs[1].axvline(1.5,color='grey',ls='--',lw=.8)
    axs[1].set(xticks=range(4),xticklabels=['GSE HSIL','GSE SCC','EBI HSIL','EBI SCC'],ylabel='Frozen program in basal-marker state',title='B  Specimen scores within one marker state',xlim=(-.35,3.45))
    axs[1].tick_params(axis='x',rotation=35)
    fig.savefig(F/'19_two_source_epithelial_states.png',dpi=300);fig.savefig(F/'19_two_source_epithelial_states.pdf');plt.close(fig)

if __name__=='__main__':main()

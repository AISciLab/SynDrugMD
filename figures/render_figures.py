from pathlib import Path
import shutil,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
O=Path(__file__).resolve().parents[1];T=O/'统计数据';F=O/'图';F.mkdir(exist_ok=True)
z=pd.read_csv(T/'analysis_records.csv.gz',keep_default_na=False,low_memory=False)
z['y']=z['Synergy?'].astype(str).str.upper().map({'TRUE':1,'FALSE':0})
primary=pd.read_csv(T/'pathway_comparisons.csv')
primary=primary[(primary.method=='cell_cluster')&(primary.family=='raw')].copy()
profile=pd.read_csv(T/'threshold_profiles_all_counts.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'savefig.facecolor':'white'})
RED='#943A2B';BLUE='#8CA0BE';GREY='#707780';NAVY='#044875'
def figsave(fig,name):
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',metadata={'Author':'','Title':name,'Creator':'Scientific analysis'})
    fig.savefig(F/(name+'.png'),dpi=230,bbox_inches='tight');plt.close(fig);shutil.copy2(F/(name+'.pdf'),O/(name+'.pdf'))
def fmtq(x):
    if not np.isfinite(x):return 'NE'
    return f'{x:.1e}' if x<.001 else f'{x:.3f}'

order=z.Tissue.value_counts().index.tolist()+['All']
fig,axes=plt.subplots(1,2,figsize=(7.1,3.8),sharey=True)
fig.subplots_adjust(left=.345,right=.91,top=.88,bottom=.14,wspace=.55)
for k,(ax,m,title) in enumerate(zip(axes,['OPC_stored','RPC_stored'],['A. OPC','B. RPC'])):
    q=primary[primary.metric==m].set_index('group').loc[order]
    finite=q.q_BH.notna();points=q.loc[finite,'difference'].to_numpy()
    lo=float(points.min());hi=float(points.max());pad=max(10,.22*(hi-lo));xmin=lo-pad;xmax=hi+pad
    if xmin>0:xmin=-10
    ax.set_xlim(xmin,xmax)
    ax.axvline(0,color='#777777',ls='--',lw=.65,zorder=0)
    for i,(_,r) in enumerate(q.iterrows()):
        if not np.isfinite(r.q_BH):
            ax.text(.5,i,'NE',ha='center',va='center',color=GREY,transform=ax.get_yaxis_transform(),fontsize=6.5)
        else:
            color=RED if r.q_BH<.05 else GREY
            left=max(xmin,r.ci_low);right=min(xmax,r.ci_high)
            ax.errorbar(r.difference,i,xerr=[[r.difference-left],[right-r.difference]],fmt='D' if r.name=='All' else 'o',markersize=3 if r.name=='All' else 2.6,capsize=1.7,lw=.8,color=color,zorder=3)
            if r.ci_low<xmin:ax.plot(xmin,i,marker='<',ms=3.5,color=color,clip_on=False)
            if r.ci_high>xmax:ax.plot(xmax,i,marker='>',ms=3.5,color=color,clip_on=False)
        ax.text(1.04,i,fmtq(r.q_BH),va='center',ha='left',fontsize=6.1,transform=ax.get_yaxis_transform())
    ax.axhline(len(order)-1.5,color='#B8BDC4',lw=.6)
    ax.set_ylim(len(order)-.4,-.6);ax.set_title(title,fontsize=8.3,loc='left',pad=9)
    ax.set_xlabel('Mean difference',fontsize=7);ax.tick_params(axis='x',labelsize=6.4);ax.tick_params(axis='y',length=0)
    ax.text(1.04,1.02,'BH q',fontsize=6.4,transform=ax.transAxes)
    if k==0:
        labels=[]
        for tissue in order:
            r=q.loc[tissue];labels.append(f'{tissue}  {int(r.n_synergy):,}/{int(r.n_non):,}')
        ax.set_yticks(np.arange(len(order)),labels,fontsize=6.3)
        ax.get_yticklabels()[-1].set_fontweight('bold')
    else:ax.tick_params(axis='y',labelleft=False)
fig.text(.015,.93,'Original label     n_syn / n_non',fontsize=6.6)
fig.text(.57,.985,'Synergistic minus non-synergistic; 95% cell-cluster intervals',ha='center',va='top',fontsize=7.2)
fig.legend(handles=[Line2D([0],[0],marker='o',color=RED,lw=0,label='BH q < 0.05'),Line2D([0],[0],marker='o',color=GREY,lw=0,label='BH q >= 0.05')],loc='lower center',bbox_to_anchor=(.57,-.005),ncol=2,fontsize=6.3,frameon=False)
figsave(fig,'Pathway_Statistics')

p=profile[(profile.metric=='RPC_stored')&(profile.threshold>=0)&profile.display_n_at_least_100].copy()
assert len(p)>0
fig,axes=plt.subplots(1,2,figsize=(7.1,2.23))
fig.subplots_adjust(left=.075,right=.985,bottom=.24,top=.85,wspace=.30)
ax=axes[0]
ax.plot(p.threshold,p.synergy_fraction,color=RED,lw=1.15,label='Synergistic')
ax.plot(p.threshold,p.non_synergy_fraction,color=BLUE,lw=1.15,label='Non-synergistic')
basefrac=float(z.y.mean());ax.axhline(basefrac,color='#555555',ls=':',lw=.85,label=f'Full-cohort synergy fraction ({100*basefrac:.2f}%)')
ax.set_ylim(0,1);ax.set_ylabel('Fraction among retained records',fontsize=7)
ax.set_title('A. Conditional class proportions',fontsize=8,loc='left')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.29),ncol=1,fontsize=5.8,frameon=False,handlelength=2)
ax=axes[1]
ax.plot(p.threshold,p.n_remaining,color=NAVY,lw=1,label='All retained')
ax.plot(p.threshold,p.n_synergy,color=RED,lw=1,label='Synergistic')
ax.plot(p.threshold,p.n_non,color=BLUE,lw=1,label='Non-synergistic')
ax.set_yscale('log');ax.set_ylabel('Retained records (log scale)',fontsize=7)
ax.set_title('B. Support at each threshold',fontsize=8,loc='left');ax.legend(loc='upper right',fontsize=6,frameon=False)
for ax in axes:
    ax.set_xlim(0,float(p.threshold.max()));ax.set_xlabel('RPC threshold (strictly greater than)',fontsize=7);ax.tick_params(labelsize=6.5);ax.grid(alpha=.12,lw=.5)
figsave(fig,'RPC_Threshold_Profiles')

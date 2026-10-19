from pathlib import Path
import sys,json,warnings,argparse
import numpy as np,pandas as pd,scipy
from scipy import stats

O=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--records',type=Path,default=O/'统计数据/analysis_records.csv.gz');ap.add_argument('--out',type=Path,default=O/'统计数据');args=ap.parse_args()
D=args.out;D.mkdir(exist_ok=True,parents=True)
z=pd.read_csv(args.records,keep_default_na=False,low_memory=False)
numeric=['OPC_stored','RPC_stored','OPC_pathway_jaccard','RPC_relation_pair_density','target_n_a','target_n_b','pathway_n_a','pathway_n_b','pathway_product','ctd_path_degree_sum_a','ctd_path_degree_sum_b']
for column in numeric:z[column]=pd.to_numeric(z[column],errors='coerce')
z['y']=z['Synergy?'].astype(str).str.upper().map({'TRUE':1,'FALSE':0})
assert z.y.notna().all()
z['source']=z.source.fillna('').replace('','(missing)')
z['Tissue']=z.Tissue.fillna('').replace('','Unknown')
z['pair_cluster']=[s+'__'+'||'.join(sorted((str(a),str(b)))) for s,a,b in zip(z.source,z['Drug1 ID'],z['Drug2 ID'])]
z['cell_cluster']=z['Cell Line name'].astype(str).str.upper()
valid_path=(z.pathway_n_a>0)&(z.pathway_n_b>0)
z.loc[~valid_path,['OPC_pathway_jaccard','RPC_relation_pair_density']]=np.nan
METRICS={'raw':['OPC_stored','RPC_stored'],'normalized':['OPC_pathway_jaccard','RPC_relation_pair_density']}
def save(d,name):d.to_csv(D/name,index=False,encoding='utf-8-sig')
def pcap(p):return max(float(p),1e-300) if np.isfinite(p) else np.nan
def bh(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p))
    if len(ix):
        order=ix[np.argsort(p[ix])];vals=p[order]*len(ix)/np.arange(1,len(ix)+1);out[order]=np.minimum.accumulate(vals[::-1])[::-1].clip(0,1)
    return out
def basic(g,metric,kind):
    g=g[np.isfinite(g[metric])];ys=g.loc[g.y==1,metric].to_numpy(float);yn=g.loc[g.y==0,metric].to_numpy(float)
    n1=len(ys);n0=len(yn);ms=ys.mean() if n1 else np.nan;mn=yn.mean() if n0 else np.nan
    out={'n_synergy':n1,'n_non':n0,'mean_synergy':ms,'mean_non':mn,'difference':ms-mn,'p':np.nan,'ci_low':np.nan,'ci_high':np.nan,'df':np.nan,'cluster_n':np.nan,'hedges_g':np.nan,'p_floored':False}
    if n1<2 or n0<2:return out
    v1=ys.var(ddof=1);v0=yn.var(ddof=1)
    pool=((n1-1)*v1+(n0-1)*v0)/(n1+n0-2)
    out['hedges_g']=(1-3/(4*(n1+n0-2)-1))*(ms-mn)/np.sqrt(pool) if pool>0 else np.nan
    if kind=='welch':
        se2=v1/n1+v0/n0
        if se2<=0:return out
        df=se2**2/((v1/n1)**2/(n1-1)+(v0/n0)**2/(n0-1));se=np.sqrt(se2)
    else:
        cluster=g['cell_cluster' if kind=='cell_cluster' else 'pair_cluster']
        ci,levels=pd.factorize(cluster,sort=True);G=len(levels)
        out['cluster_n']=G
        if G<2:return out
        y=g[metric].to_numpy(float);s=g.y.to_numpy(int)
        influence=np.where(s==1,(y-ms)/n1,-(y-mn)/n0)
        u=np.bincount(ci,weights=influence,minlength=G)
        se2=G/(G-1)*(len(g)-1)/(len(g)-2)*np.sum(u*u)
        if se2<=0:return out
        df=G-1;se=np.sqrt(se2)
    p=2*stats.t.sf(abs((ms-mn)/se),df);crit=stats.t.ppf(.975,df)
    out.update({'p':pcap(p),'p_floored':bool(p<1e-300),'ci_low':ms-mn-crit*se,'ci_high':ms-mn+crit*se,'df':df})
    return out

records=[];groups=[('All',z)]+list(z.groupby('Tissue',sort=True))
for family,metrics in METRICS.items():
    for kind in ['welch','cell_cluster','source_pair_cluster']:
        batch=[]
        for tissue,g in groups:
            for m in metrics:
                a=basic(g,m,kind);a.update({'family':family,'method':kind,'group':tissue,'metric':m});batch.append(a)
        qq=bh([a['p'] for a in batch]);estimable=int(np.isfinite([a['p'] for a in batch]).sum())
        for a,q in zip(batch,qq):a.update({'q_BH':q,'family_tests_planned':len(batch),'family_tests_estimable':estimable})
        records.extend(batch)
tt=pd.DataFrame(records);save(tt,'pathway_comparisons.csv')
save(tt[tt.group=='All'],'overall_effects.csv')
sig=tt[tt.group!='All'].groupby(['family','method','metric']).agg(estimable_tissue_labels=('q_BH',lambda s:int(s.notna().sum())),q_below_005=('q_BH',lambda s:int((s<.05).sum()))).reset_index()
save(sig,'significance_summary.csv')

# Fixed-sample check: common affine transformations leave the t statistic invariant.
x=z.loc[z.y==1,'OPC_stored'].to_numpy(float);y=z.loc[z.y==0,'OPC_stored'].to_numpy(float)
allv=z.OPC_stored.to_numpy(float);low=allv.min();span=allv.max()-low;mean=allv.mean();std=allv.std(ddof=1)
t0=stats.ttest_ind(x,y,equal_var=True);t1=stats.ttest_ind((x-low)/span,(y-low)/span,equal_var=True);t2=stats.ttest_ind((x-mean)/std,(y-mean)/std,equal_var=True)
affine={'raw_t':float(t0.statistic),'minmax_t':float(t1.statistic),'zscore_t':float(t2.statistic),'max_abs_t_difference':float(max(abs(t0.statistic-t1.statistic),abs(t0.statistic-t2.statistic))),
 'interpretation':'Common Min-Max and global Z-score are affine rescalings. On identical observations they do not adjust target/pathway-size confounding.'}
(D/'affine_invariance.json').write_text(json.dumps(affine,ensure_ascii=False,indent=2),encoding='utf8')

def targetbin(v):return '1' if v==1 else '2-5' if v<=5 else '6+'
def log4bin(v):return int(np.floor(np.log(max(float(v),1))/np.log(4)))
comp=z[valid_path].copy()
comp['target_bin']=[targetbin(min(a,b))+'/'+targetbin(max(a,b)) for a,b in zip(comp.target_n_a,comp.target_n_b)]
comp['path_bin']=comp.pathway_product.map(log4bin)
comp['degree_bin']=(comp.ctd_path_degree_sum_a*comp.ctd_path_degree_sum_b).map(log4bin)
schemes={'source':['source'],'source_targets':['source','target_bin'],'source_targets_pathways':['source','target_bin','path_bin'],'source_targets_pathways_degree':['source','target_bin','path_bin','degree_bin']}
mm=['OPC_stored','RPC_stored','OPC_pathway_jaccard','RPC_relation_pair_density']
adjusted=[];strata_rows=[]
for scheme,fields in schemes.items():
    q=comp.copy();q['stratum']=q[fields].astype(str).agg('|'.join,axis=1)
    for m in mm:q[m+'_sum']=q[m]
    agg=q.groupby(['pair_cluster','source','stratum','y'],sort=True).agg(n=('y','size'),**{m:(m+'_sum','sum') for m in mm}).reset_index()
    ci,cn=pd.factorize(agg.pair_cluster,sort=True);si,sn=pd.factorize(agg.stratum,sort=True)
    ns=len(sn);nc=len(cn);label=agg.y.to_numpy(int);nn=agg.n.to_numpy(float)
    sy=np.bincount(si,weights=nn*(label==1),minlength=ns);no=np.bincount(si,weights=nn*(label==0),minlength=ns)
    csy=np.bincount(si,weights=(label==1),minlength=ns);cno=np.bincount(si,weights=(label==0),minlength=ns)
    support=(sy>0)&(no>0)&(csy>=5)&(cno>=5)
    weights=(sy+no)*support;weights/=weights.sum()
    source_c=agg[['pair_cluster','source']].drop_duplicates().set_index('pair_cluster').loc[cn,'source'].to_numpy()
    cgroups=[np.flatnonzero(source_c==s) for s in sorted(set(source_c))]
    vals={m:[] for m in mm};points={}
    for m in mm:
        sums=agg[m].to_numpy(float);ss=np.bincount(si,weights=sums*(label==1),minlength=ns);nsum=np.bincount(si,weights=sums*(label==0),minlength=ns)
        ps=float(np.sum(weights[support]*ss[support]/sy[support]));pn=float(np.sum(weights[support]*nsum[support]/no[support]));points[m]=(ps,pn)
        for i in np.flatnonzero(support):strata_rows.append({'scheme':scheme,'stratum':sn[i],'metric':m,'n_synergy':sy[i],'n_non':no[i],'source_pair_clusters_synergy':csy[i],'source_pair_clusters_non':cno[i],'weight':weights[i],'mean_synergy':ss[i]/sy[i],'mean_non':nsum[i]/no[i]})
    rng=np.random.default_rng(20261019)
    for _ in range(500):
        ix=np.concatenate([rng.choice(g,len(g),replace=True) for g in cgroups]);mul=np.bincount(ix,minlength=nc)[ci]
        bs=np.bincount(si,weights=mul*nn*(label==1),minlength=ns);bn=np.bincount(si,weights=mul*nn*(label==0),minlength=ns)
        if np.any((bs[support]==0)|(bn[support]==0)):continue
        for m in mm:
            sums=agg[m].to_numpy(float);ss=np.bincount(si,weights=mul*sums*(label==1),minlength=ns);nsum=np.bincount(si,weights=mul*sums*(label==0),minlength=ns)
            vals[m].append(float(np.sum(weights[support]*(ss[support]/bs[support]-nsum[support]/bn[support]))))
    for m in mm:
        lo,hi=np.quantile(vals[m],[.025,.975]) if vals[m] else (np.nan,np.nan)
        adjusted.append({'scheme':scheme,'metric':m,'eligible_records':len(comp),'common_support_records':int(((sy+no)*support).sum()),'common_support_fraction':float(((sy+no)*support).sum()/len(comp)),
          'supported_strata':int(support.sum()),'source_pair_cluster_n':nc,'standardized_mean_synergy':points[m][0],'standardized_mean_non':points[m][1],'standardized_difference':points[m][0]-points[m][1],
          'bootstrap_95_low':lo,'bootstrap_95_high':hi,'bootstrap_requested':500,'bootstrap_valid':len(vals[m]),'seed':20261019})
    print('Standardization completed',scheme,'common support',float(((sy+no)*support).sum()/len(comp)),flush=True)
ad=pd.DataFrame(adjusted);save(ad,'size_standardization.csv');save(pd.DataFrame(strata_rows),'stratum_support.csv')

eq=[]
for c,g in z.groupby('cell_cluster'):
    for m in mm:
        a=g.loc[g.y==1,m].dropna();b=g.loc[g.y==0,m].dropna()
        if len(a) and len(b):eq.append({'cell':c,'metric':m,'difference':a.mean()-b.mean(),'n_synergy':len(a),'n_non':len(b)})
eqd=pd.DataFrame(eq);save(eqd,'within_cell_differences.csv')
eqsum=eqd.groupby('metric').agg(equal_cell_difference=('difference','mean'),cells_with_both_labels=('cell','size')).reset_index();save(eqsum,'equal_cell_differences.csv')

summary={'software':{'scipy':scipy.__version__,'pandas':pd.__version__,'numpy':np.__version__},'cohort_rows':len(z),'target_pathway_complete_case_rows':len(comp),
 'n_synergy':int((z.y==1).sum()),'n_non':int((z.y==0).sum()),'cell_clusters':z.cell_cluster.nunique(),'source_pair_clusters':z.pair_cluster.nunique(),
 'metrics':METRICS,'overall':tt[tt.group=='All'].replace({np.nan:None}).to_dict('records'),'stratified':ad.to_dict('records'),'equal_cell':eqsum.to_dict('records'),
 'affine_invariance':affine,'inference_limits':['Normalized metrics are newly constructed sensitivity analyses, not results proven to have been performed in the original notebooks.','Relation-density normalization uses the corrected direct-pathway relation count.','One-way clustering and coarse common-support standardization do not remove every source, drug, cell, annotation, or selection confound.','No degree-preserving network randomization, independent prediction validation, or causal experiment is claimed.','BH families are separate by test method and metric family; exploratory method alternatives do not constitute independent confirmations.']}
(D/'reproduced_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'rows':len(z),'complete':len(comp),'overall':tt[(tt.group=='All')&(tt.method=='cell_cluster')].replace({np.nan:None}).to_dict('records'),'adjusted':ad.to_dict('records'),'affine':affine},ensure_ascii=False,indent=2),flush=True)

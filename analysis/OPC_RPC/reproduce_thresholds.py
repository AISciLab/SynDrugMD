from pathlib import Path
import numpy as np,pandas as pd
O=Path(__file__).resolve().parents[1];T=O/'统计数据'
z=pd.read_csv(T/'analysis_records.csv.gz',keep_default_na=False,low_memory=False)
z['y']=z['Synergy?'].astype(str).str.upper().map({'TRUE':1,'FALSE':0})
def save(df,name):df.to_csv(T/name,index=False,encoding='utf-8-sig')
# Exact strict-threshold counts, including a -1 row representing all nonnegative scores.
profiles=[]
for metric in ['OPC_stored','RPC_stored']:
    x=z[metric].to_numpy(int);syn=np.sort(x[z.y.to_numpy()==1]);non=np.sort(x[z.y.to_numpy()==0])
    for tau in range(-1,int(x.max())+1):
        ns=len(syn)-int(np.searchsorted(syn,tau,side='right'));nn=len(non)-int(np.searchsorted(non,tau,side='right'));n=ns+nn
        profiles.append({'metric':metric,'threshold':tau,'strict_comparison':'>','n_synergy':ns,'n_non':nn,'n_remaining':n,'synergy_fraction':ns/n if n else np.nan,'non_synergy_fraction':nn/n if n else np.nan,'display_n_at_least_100':n>=100})
profile=pd.DataFrame(profiles);save(profile,'threshold_profiles_all_counts.csv')
for metric in ['OPC_stored','RPC_stored']:
    p=profile[profile.metric==metric]
    assert (p.n_remaining.diff().dropna()<=0).all()
    f=p[p.n_remaining>0]
    assert np.allclose(f.synergy_fraction+f.non_synergy_fraction,1)
    for tau in [-1,0,20,100,500,1000,int(z[metric].max())]:
        if tau>z[metric].max():continue
        r=p[p.threshold==tau].iloc[0];keep=z[metric]>tau
        assert int(keep.sum())==r.n_remaining and int(z.loc[keep,'y'].sum())==r.n_synergy


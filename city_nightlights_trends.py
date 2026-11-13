import sys
DIR = sys.argv[1] if len(sys.argv) > 1 else "."   # folder with the CSV files
import pandas as pd, numpy as np, glob, json
from scipy.stats import norm, rankdata, theilslopes
def hr_mk(t,y,maxlag=None):
    y=np.asarray(y,float); n=len(y)
    sl=theilslopes(y,t)[0]
    det=y-sl*t; r=rankdata(det)
    S=sum(np.sign(y[j]-y[i]) for i in range(n-1) for j in range(i+1,n))
    v=n*(n-1)*(2*n+5)/18
    maxlag=maxlag or n-1; rr=r-r.mean(); den=(rr**2).sum(); a=[]
    lim=norm.ppf(0.975)/np.sqrt(n)
    f=0
    for k in range(1,min(maxlag,n-2)+1):
        rho=(rr[:-k]*rr[k:]).sum()/den
        if abs(rho)>lim: f+=(n-k)*(n-k-1)*(n-k-2)*rho
    cf=max(1.0,1+2/(n*(n-1)*(n-2))*f)
    z=(S-np.sign(S))/np.sqrt(v*cf) if S else 0
    return sl, 2*(1-norm.cdf(abs(z))), cf
def bh(p):
    p=np.asarray(p); o=np.argsort(p); q=np.empty_like(p); m=len(p); prev=1
    for rank,i in zip(range(m,0,-1),o[::-1]): prev=min(prev,p[i]*m/rank); q[i]=prev
    return q
D={}
for f in sorted(glob.glob(DIR+'/'+'*nightlights_*.csv')):
    n=f.split('nightlights_')[1][:-4]; d=pd.read_csv(f); d['date']=pd.to_datetime(d['date']); D[n]=d.set_index('date')['mean']
df=pd.DataFrame(D)
def run(sub,label):
    rows=[]
    for c in sub.columns:
        s=sub[c]; clim=s.groupby(s.index.month).transform('mean'); an=s-clim
        t=(s.index.year+(s.index.month-1)/12).values; t=t-t.mean()
        sl,p,cf=hr_mk(t,an.values,24)
        rows.append(dict(city=c,mean=s.mean(),slope_dec=sl*10,pct_dec=sl*10/s.mean()*100,p=p,infl=cf))
    r=pd.DataFrame(rows); r['q']=bh(r.p.values); return r
full=run(df,'full'); print(full.round(4).to_string())
post=run(df[df.index>='2014-01'],'2014+'); print('\n2014 onward'); print(post.round(4).to_string())
no20=run(df[(df.index.year!=2020)],'no2020'); print('\nno 2020'); print(no20.round(4).to_string())
ann=df.groupby(df.index.year).mean(); ann=ann[ann.index<=2025]  # full years only
print('\nlast partial year months:',df[df.index.year==2026].shape[0])
json.dump({'years':list(map(int,ann.index)),'ann':{c:[round(v,1) for v in ann[c]] for c in ann.columns},
           'full':full.round(4).to_dict('records'),'post':post.round(4).to_dict('records')},open('nl.json','w'),ensure_ascii=False)

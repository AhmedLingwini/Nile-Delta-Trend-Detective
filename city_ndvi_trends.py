import sys
DIR = sys.argv[1] if len(sys.argv) > 1 else "."   # folder with the CSV files
import pandas as pd, numpy as np, glob, json
from scipy.stats import norm, rankdata, theilslopes, kendalltau
def hr_mk(t,y,maxlag=24):
    y=np.asarray(y,float); n=len(y); sl=theilslopes(y,t)[0]
    r=rankdata(y-sl*t); S=np.sign(y[None,:]-y[:,None])[np.triu_indices(n,1)].sum()
    v=n*(n-1)*(2*n+5)/18; rr=r-r.mean(); den=(rr**2).sum(); lim=norm.ppf(.975)/np.sqrt(n); f=0
    for k in range(1,min(maxlag,n-2)+1):
        rho=(rr[:-k]*rr[k:]).sum()/den
        if abs(rho)>lim: f+=(n-k)*(n-k-1)*(n-k-2)*rho
    cf=max(1.0,1+2/(n*(n-1)*(n-2))*f); z=(S-np.sign(S))/np.sqrt(v*cf)
    return sl,2*(1-norm.cdf(abs(z))),cf
def bh(p):
    p=np.asarray(p); o=np.argsort(p); q=np.empty_like(p); m=len(p); prev=1
    for rank,i in zip(range(m,0,-1),o[::-1]): prev=min(prev,p[i]*m/rank); q[i]=prev
    return q
nd=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(DIR+'/'+'ndvi_20*.csv'))])
nd['dt']=pd.to_datetime(nd['date']); nd['year']=nd.dt.dt.year; nd['month']=nd.dt.dt.month
reg=nd.groupby('year')['mean'].mean(); print('regional',reg.round(3).iloc[[0,-1]].to_dict(), 'rows',len(nd))
C={'Alexandria':(29.92,31.20),'Mansoura':(31.38,31.04),'Damietta':(31.82,31.42),'Cairo':(31.24,30.04),'Ismailia':(32.27,30.59)}
out={'reg':[round(v,4) for v in reg.values],'years':list(map(int,reg.index)),'city':{}}
rows=[]
for n,(lo,la) in C.items():
    km=np.hypot((nd.lon-lo)*np.cos(np.radians(la))*111,(nd.lat-la)*111)
    s=nd[km<=15]; cells=s.cell_id.nunique()
    ann=s.groupby('year')['mean'].mean()
    m=s.groupby('dt')['mean'].mean(); an=m-m.groupby(m.index.month).transform('mean')
    t=(m.index.year+(m.index.month-1)/12).values; t=t-t.mean()
    sl,p,cf=hr_mk(t,an.values)
    sa=theilslopes(ann.values,ann.index.values)[0]*10
    rows.append((n,cells,sl*10,p)); out['city'][n]={'ndvi':[round(v,4) for v in ann.values],'cells':int(cells),'ndvi_dec':round(sl*10,4),'ndvi_p':float(p),'ndvi_mean':round(float(ann.mean()),3)}
q=bh([r[3] for r in rows])
for (n,c,sl,p),qq in zip(rows,q): out['city'][n]['ndvi_q']=float(qq); print(n,'cells',c,'slope/dec %.4f p %.4f q %.4f'%(sl,p,qq))
# regional stats like script
yy=reg.index.values; tau,p=kendalltau(yy,reg.values); print('regional MK p',p,'Sen/dec',theilslopes(reg.values,yy)[0]*10)
out['reg_stats']={'sen_dec':round(theilslopes(reg.values,yy)[0]*10,4),'p':float(p)}
json.dump(out,open('ndvi_out.json','w'))

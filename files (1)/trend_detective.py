"""
Earth System Trend Detective - Nile Delta NDVI (2001-2025) + NASA POWER climate
Usage:  python trend_detective.py --ndvi_dir delta_trend_detective --power_dir . --out out
Inputs : ndvi_*.csv  (cell_id,lon,lat,date,mean)   and   POWER_Point_Monthly_*.csv
Outputs: figures/*.png, cell_trends.csv, summary.txt
"""
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--ndvi_dir", default="delta_trend_detective")
ap.add_argument("--power_dir", default=".")
ap.add_argument("--out", default="out")
a = ap.parse_args()
os.makedirs(f"{a.out}/figures", exist_ok=True)
log = []
def say(s): print(s); log.append(s)

# ---------- 1. Load NDVI ----------
nd = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{a.ndvi_dir}/ndvi_*.csv"))])
nd["date"] = pd.to_datetime(nd["date"]); nd["year"] = nd.date.dt.year; nd["month"] = nd.date.dt.month
cells = nd.groupby("cell_id")["mean"].apply(lambda s: s.notna().mean())
valid = cells[cells >= 0.5].index           # drop all-NaN cells (sea / lakes)
say(f"Cells total {nd.cell_id.nunique()}, usable (>=50% data) {len(valid)}")
nd = nd[nd.cell_id.isin(valid)]

# ---------- 2. Regional trend ----------
ann = nd.groupby("year")["mean"].mean()
lr = stats.linregress(ann.index, ann.values)
say(f"Regional mean NDVI: {ann.iloc[0]:.3f} (2001) -> {ann.iloc[-1]:.3f} (2025); "
    f"slope {lr.slope:+.5f}/yr, R2={lr.rvalue**2:.2f}, p={lr.pvalue:.1e}")
mk = stats.kendalltau(ann.index, ann.values); say(f"Mann-Kendall tau={mk[0]:.2f}, p={mk[1]:.1e}")

# ---------- 3. Per-cell trends (Theil-Sen + Mann-Kendall on annual means) ----------
pa = nd.groupby(["cell_id", "lon", "lat", "year"])["mean"].mean().reset_index()
rows = []
for (c, lo, la), g in pa.groupby(["cell_id", "lon", "lat"]):
    g = g.dropna()
    if len(g) < 15: continue
    tau, p = stats.kendalltau(g.year, g["mean"])
    rows.append((c, lo, la, stats.theilslopes(g["mean"], g.year)[0], tau, p, g["mean"].mean()))
tr = pd.DataFrame(rows, columns=["cell_id", "lon", "lat", "slope_per_yr", "tau", "p", "mean_ndvi"])
tr["class"] = np.where(tr.p >= .05, "no significant trend", np.where(tr.slope_per_yr > 0, "greening", "browning"))
tr.to_csv(f"{a.out}/cell_trends.csv", index=False)
say("Per-cell classes:\n" + tr["class"].value_counts().to_string())
say("Top browning cells:\n" + tr.nsmallest(5, "slope_per_yr")[["cell_id","lon","lat","slope_per_yr","p"]].round(4).to_string(index=False))
say("Top greening cells:\n" + tr.nlargest(5, "slope_per_yr")[["cell_id","lon","lat","slope_per_yr","p"]].round(4).to_string(index=False))


# ---------- 3b. CONFIDENCE METER: "is the trend real?" ----------
# (a) monthly de-seasonalised anomalies  (b) Mann-Kendall with Hamed-Rao autocorrelation correction
# (c) Benjamini-Hochberg FDR across all cells  (d) agreement with the independent annual-mean test
def mk_hamed_rao(x):
    x = np.asarray(x, float); t = np.arange(len(x)); m = ~np.isnan(x); x, t = x[m], t[m]; n = len(x)
    if n < 60: return None
    S = np.sum(np.sign(x[None, :] - x[:, None])[np.triu_indices(n, 1)])
    sl = stats.theilslopes(x, t)[0]; rk = stats.rankdata(x - sl * t); rk = rk - rk.mean()
    ac = np.array([np.sum(rk[:-k] * rk[k:]) / np.sum(rk * rk) for k in range(1, n - 1)])
    ks = np.arange(1, n - 1); sg = np.abs(ac) > 1.96 / np.sqrt(n)
    vif = max(1 + 2 / (n*(n-1)*(n-2)) * np.sum((n-ks[sg])*(n-ks[sg]-1)*(n-ks[sg]-2)*ac[sg]), 1.0)
    v0 = n*(n-1)*(2*n+5)/18; z = (S - np.sign(S)) / np.sqrt(v0 * vif)
    pn = 2*(1 - stats.norm.cdf(abs((S - np.sign(S)) / np.sqrt(v0))))
    return sl * 12, 2*(1 - stats.norm.cdf(abs(z))), vif, ac[0], pn
def bh_fdr(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); r = np.minimum.accumulate((p[o]*n/np.arange(1, n+1))[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(r, 1); return q
vlim = np.nanpercentile(np.abs(tr.slope_per_yr), 98)
Pm = nd.pivot(index="date", columns="cell_id", values="mean").sort_index()
An = Pm - Pm.groupby(Pm.index.month).transform("mean")
cr = []
for c in An.columns:
    r = mk_hamed_rao(An[c].values)
    if r: cr.append((c,) + r)
cr = pd.DataFrame(cr, columns=["cell_id", "slope_monthly_per_yr", "p_autocorr", "autocorr_vif", "lag1", "p_naive_monthly"])
cr["q_fdr"] = bh_fdr(cr.p_autocorr)
tr = tr.merge(cr, on="cell_id", how="left")
agree = (np.sign(tr.slope_per_yr) == np.sign(tr.slope_monthly_per_yr)) & (tr.p < .05)
tr["confidence"] = np.where((tr.q_fdr < .05) & agree, "High", np.where(tr.q_fdr < .05, "Medium", "Low"))
tr["direction"] = np.where(tr.slope_per_yr > 0, "greening", "browning")
tr.to_csv(f"{a.out}/cell_trends.csv", index=False)
say(f"Autocorrelation: median lag-1 = {tr.lag1.median():.2f}, median variance inflation = {tr.autocorr_vif.median():.1f}x")
say(f"Significant cells (p<0.05): naive monthly MK {(tr.p_naive_monthly<.05).sum()} | autocorr-corrected {(tr.p_autocorr<.05).sum()} | + FDR {(tr.q_fdr<.05).sum()} of {len(tr)}")
say("Confidence x direction:\n" + pd.crosstab(tr.confidence, tr.direction).to_string())
cp = tr.pivot(index="lat", columns="lon", values="confidence").replace({"High": 2, "Medium": 1, "Low": 0}).astype(float)
dp = tr.pivot(index="lat", columns="lon", values="slope_per_yr")
fig, ax = plt.subplots(figsize=(7.5, 5.5)); show = dp.where(cp == 2)
ax.imshow(np.where(cp.values < 2, 0.5, np.nan), origin="lower", cmap="Greys", vmin=0, vmax=1, extent=[dp.columns.min()-.05, dp.columns.max()+.05, dp.index.min()-.05, dp.index.max()+.05])
im = ax.imshow(show.values, origin="lower", cmap="BrBG", vmin=-vlim, vmax=vlim, extent=[dp.columns.min()-.05, dp.columns.max()+.05, dp.index.min()-.05, dp.index.max()+.05])
ax.set(title="High-confidence trends only (grey = medium/low confidence)", xlabel="Longitude", ylabel="Latitude"); plt.colorbar(im, label="NDVI change per year")
fig.tight_layout(); fig.savefig(f"{a.out}/figures/5_confidence_map.png"); plt.close()

# ---------- 4. Climate (NASA POWER) ----------
cl = []
for f in glob.glob(f"{a.power_dir}/*POWER_Point_Monthly*.csv"):
    d = pd.read_csv(f, skiprows=10); d = d[d.YEAR.between(2001, 2025)]
    cl.append(d.assign(src=os.path.basename(f)))
cl = pd.concat(cl)
t2m = cl[cl.PARAMETER == "T2M"].groupby("YEAR")["ANN"].mean()
pr = cl[cl.PARAMETER == "PRECTOTCORR"].groupby("YEAR")["ANN"].mean()
for name, s in [("T2M (C)", t2m), ("Precip (mm/day)", pr)]:
    l = stats.linregress(s.index, s.values)
    say(f"{name}: slope {l.slope:+.4f}/yr  -> {l.slope*24:+.2f} over 2001-2025, p={l.pvalue:.1e}")
# correlation: raw (both trend -> spurious risk) vs de-trended
def detr(s): return pd.Series(s.values - np.polyval(np.polyfit(s.index, s.values, 1), s.index), index=s.index)
for name, s in [("T2M", t2m), ("Precip", pr)]:
    r_raw = stats.pearsonr(ann, s.reindex(ann.index))
    r_det = stats.pearsonr(detr(ann), detr(s.reindex(ann.index)))
    say(f"NDVI vs {name}: raw r={r_raw[0]:+.2f} (p={r_raw[1]:.3f}) | de-trended r={r_det[0]:+.2f} (p={r_det[1]:.3f})")

# ---------- 5. Figures ----------
plt.rcParams.update({"figure.dpi": 130, "axes.spines.top": False, "axes.spines.right": False})
# 5a regional + climate
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(ann.index, ann.values, "o-", c="#2e7d32"); ax[0].plot(ann.index, lr.intercept + lr.slope*ann.index, "--", c="k")
ax[0].set(title=f"Nile Delta mean NDVI ({lr.slope*10:+.3f} per decade)", xlabel="Year", ylabel="NDVI")
ax[1].plot(t2m.index, t2m.values, "o-", c="#d84315", label="T2M (°C)"); ax[1].set(title="Climate at 5 POWER points", xlabel="Year", ylabel="Temperature (°C)")
b = ax[1].twinx(); b.plot(pr.index, pr.values, "s-", c="#1565c0", label="Precip"); b.set_ylabel("Precip (mm/day)")
fig.tight_layout(); fig.savefig(f"{a.out}/figures/1_regional_trend.png"); plt.close()
# 5b trend map
piv = tr.pivot(index="lat", columns="lon", values="slope_per_yr"); sig = tr.pivot(index="lat", columns="lon", values="p") < .05
fig, ax = plt.subplots(figsize=(7.5, 5.5)); v = np.nanpercentile(np.abs(tr.slope_per_yr), 98)
im = ax.imshow(piv.values, origin="lower", cmap="BrBG", vmin=-v, vmax=v, extent=[piv.columns.min()-.05, piv.columns.max()+.05, piv.index.min()-.05, piv.index.max()+.05])
yy, xx = np.where(~sig.values & ~np.isnan(piv.values)); ax.scatter(piv.columns.values[xx], piv.index.values[yy], s=3, c="k", alpha=.35, label="not significant")
ax.set(title="NDVI trend per cell (Theil–Sen, 2001–2025)", xlabel="Longitude", ylabel="Latitude"); ax.legend(loc="lower left", fontsize=7)
plt.colorbar(im, label="NDVI change per year"); fig.tight_layout(); fig.savefig(f"{a.out}/figures/2_trend_map.png"); plt.close()
# 5c seasonal cycle early vs late
seas = nd.assign(p=np.where(nd.year <= 2005, "2001–2005", np.where(nd.year >= 2021, "2021–2025", None))).dropna(subset=["p"]).groupby(["p", "month"])["mean"].mean().unstack(0)
fig, ax = plt.subplots(figsize=(6.5, 4)); seas.plot(ax=ax, marker="o"); ax.set(title="Seasonal NDVI cycle: early vs recent", xlabel="Month", ylabel="NDVI")
fig.tight_layout(); fig.savefig(f"{a.out}/figures/3_seasonal.png"); plt.close()
# 5d class counts
fig, ax = plt.subplots(figsize=(5, 3.5)); tr["class"].value_counts().reindex(["greening","no significant trend","browning"]).plot.bar(ax=ax, color=["#2e7d32","#9e9e9e","#8d6e63"])
ax.set(title="Cells by NDVI trend class", ylabel="cells"); plt.xticks(rotation=15); fig.tight_layout(); fig.savefig(f"{a.out}/figures/4_classes.png"); plt.close()

open(f"{a.out}/summary.txt", "w").write("\n".join(log))

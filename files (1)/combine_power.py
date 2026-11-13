#!/usr/bin/env python3
"""Validate, combine and analyse NASA POWER monthly point CSVs (the ones downloaded from the POWER Data Access Viewer).

Usage:
    python combine_power.py --data PATH_TO_FOLDER_WITH_CSVS --out results_power

Does three things:
  1. Validates every file (header, 2001-2025 complete, no -999, ANN column matches the months,
     plausible ranges, seasonality sanity, no two files carrying the same series).
  2. Combines them into one tidy CSV (power_combined_long.csv).
  3. Runs the same trend method as the Delta pipeline (Mann-Kendall + Sen's slope, autocorrelation-adjusted p,
     Benjamini-Hochberg FDR over all tests of the same variable) and writes an Excel report.

POWER temperature/precipitation come from MERRA-2 (a reanalysis MODEL, ~0.5 x 0.625 deg grid), not from
weather stations. A point is the average of a ~50 km cell.
"""
import argparse
import calendar
import glob
import io
import os
import re

import numpy as np
import pandas as pd

from config import SEASONS
from trend_analysis import bh_qvalues, mann_kendall, season_series

MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
PLACES = {  # (lat, lon): (arabic name, english name)
    (31.20, 29.90): ("الإسكندرية / غرب الدلتا", "Alexandria / West Delta"),
    (31.04, 31.38): ("المنصورة / قلب الدلتا", "Mansoura / Delta core"),
    (31.42, 31.82): ("دمياط / شرق الدلتا", "Damietta / East Delta"),
    (30.04, 31.24): ("القاهرة", "Cairo"),
    (30.59, 32.27): ("الإسماعيلية", "Ismailia"),
}
UNITS = {"T2M": "°C", "PRECTOTCORR": "mm/day"}
PNAME = {"T2M": "حرارة الهواء 2م", "PRECTOTCORR": "الأمطار (مصحّحة)"}
SEASON_AR = {"annual": "سنوي", "DJF": "شتاء (ديسمبر–فبراير)", "MAM": "ربيع (مارس–مايو)",
             "JJA": "صيف (يونيو–أغسطس)", "SON": "خريف (سبتمبر–نوفمبر)"}


def place_of(lat, lon):
    for (la, lo), (ar, en) in PLACES.items():
        if abs(la - lat) < 0.05 and abs(lo - lon) < 0.05:
            return ar, en
    return f"{lat:.2f}N {lon:.2f}E", f"{lat:.2f}N {lon:.2f}E"


def parse_power_csv(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    end = next(i for i, l in enumerate(lines) if l.strip() == "-END HEADER-")
    head = lines[:end]
    joined = " ".join(head)
    m = re.search(r"Latitude\s+(-?[\d.]+)\s+Longitude\s+(-?[\d.]+)", joined)
    if not m:
        raise ValueError("no Latitude/Longitude in header")
    lat, lon = float(m.group(1)), float(m.group(2))
    e = re.search(r"=\s*([\d.]+)\s*meters", joined)
    elev = float(e.group(1)) if e else np.nan
    miss = re.search(r"missing source data.*?:\s*(-?\d+)", joined)
    fill = float(miss.group(1)) if miss else -999.0
    df = pd.read_csv(io.StringIO("\n".join(lines[end + 1:])))
    return dict(lat=lat, lon=lon, elev=elev, fill=fill, df=df, path=path)


def validate(rec):
    df, fill = rec["df"], rec["fill"]
    rows, ok = [], True

    def add(check, passed, detail=""):
        nonlocal ok
        ok &= bool(passed)
        rows.append(dict(الفحص=check, النتيجة="✔" if passed else "✘", تفاصيل=str(detail)))

    params = sorted(df["PARAMETER"].unique())
    add("المتغيرات موجودة (T2M + PRECTOTCORR)", {"T2M", "PRECTOTCORR"} <= set(params), ", ".join(params))
    for p in params:
        d = df[df["PARAMETER"] == p]
        yrs = sorted(d["YEAR"].tolist())
        add(f"{p}: السنين 2001–2025 كاملة وبدون تكرار", yrs == list(range(2001, 2026)),
            f"{len(yrs)} سنة، من {yrs[0]} إلى {yrs[-1]}")
        vals = d[MONTHS].to_numpy(float)
        n_fill = int((vals == fill).sum())
        add(f"{p}: مفيش قيم مفقودة ({int(fill)})", n_fill == 0 and np.isfinite(vals).all(), f"{n_fill} قيمة مفقودة")
        v = np.where(vals == fill, np.nan, vals)
        # POWER's ANN is the day-weighted mean of the months (not the plain mean of 12 numbers)
        days = np.array([[calendar.monthrange(int(y), m)[1] for m in range(1, 13)] for y in d["YEAR"]], float)
        recomputed = (v * days).sum(axis=1) / days.sum(axis=1)
        diff = np.abs(recomputed - d["ANN"].to_numpy(float))
        add(f"{p}: عمود ANN = متوسط الشهور الموزون بعدد الأيام", diff.max() <= 0.01,
            f"أقصى فرق {diff.max():.3f} (التقريب في الملف 0.01)")
        if p == "T2M":
            add("T2M: القيم الشهرية في مدى منطقي (5–36 °C)", np.nanmin(v) > 5 and np.nanmax(v) < 36,
                f"من {np.nanmin(v):.1f} إلى {np.nanmax(v):.1f}")
            clim = np.nanmean(v, axis=0)
            add("T2M: الأسخن يوليو/أغسطس والأبرد يناير/فبراير",
                int(np.argmax(clim)) in (6, 7) and int(np.argmin(clim)) in (0, 1),
                f"الأسخن {MONTHS[int(np.argmax(clim))]}، الأبرد {MONTHS[int(np.argmin(clim))]}")
        else:
            add("PRECTOTCORR: غير سالبة وأقل من 15 mm/day", np.nanmin(v) >= 0 and np.nanmax(v) < 15,
                f"من {np.nanmin(v):.2f} إلى {np.nanmax(v):.2f}")
            dry = np.nanmean(v[:, 5:8])
            add("PRECTOTCORR: الصيف شبه جاف (يونيو–أغسطس)", dry < 0.1, f"متوسط الصيف {dry:.3f} mm/day")
    return rows, ok


def series_array(df, param):
    d = df[df["PARAMETER"] == param].sort_values("YEAR")
    arr = d[MONTHS].to_numpy(float)
    return d["YEAR"].to_numpy(int), arr


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="folder with POWER_Point_Monthly_*.csv files")
    ap.add_argument("--out", default="results_power")
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.data, "*POWER_Point_Monthly*.csv")))
    if not files:
        raise SystemExit("no POWER_Point_Monthly*.csv files found")
    recs = [parse_power_csv(f) for f in files]

    # ---- 1. validation
    val_rows, all_ok = [], True
    for r in recs:
        ar, en = place_of(r["lat"], r["lon"])
        r["ar"], r["en"] = ar, en
        rows, ok = validate(r)
        all_ok &= ok
        for x in rows:
            val_rows.append(dict(النقطة=ar, **x))
    # no two files with the same series (two close points can fall in the same ~50 km grid cell)
    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            a = series_array(recs[i]["df"], "T2M")[1]
            b = series_array(recs[j]["df"], "T2M")[1]
            same = np.allclose(a, b)
            all_ok &= not same
            val_rows.append(dict(النقطة=f"{recs[i]['ar']} ↔ {recs[j]['ar']}", الفحص="سلسلتين مختلفتين (مش نفس الخلية)",
                                 النتيجة="✘" if same else "✔",
                                 تفاصيل=f"أقصى فرق T2M {np.abs(a - b).max():.2f} °C"))
    validation = pd.DataFrame(val_rows)
    print(validation.to_string(index=False))
    print("\nALL FILES VALID" if all_ok else "\nSOME CHECKS FAILED - see above")

    # ---- 2. combine
    long_parts = []
    for r in recs:
        for p in ("T2M", "PRECTOTCORR"):
            yrs, arr = series_array(r["df"], p)
            arr = np.where(arr == r["fill"], np.nan, arr)
            for mi, mname in enumerate(MONTHS):
                long_parts.append(pd.DataFrame(dict(
                    point_ar=r["ar"], point_en=r["en"], lat=r["lat"], lon=r["lon"], elevation_m=r["elev"],
                    parameter=p, unit=UNITS[p], year=yrs, month=mi + 1, value=arr[:, mi])))
    long = pd.concat(long_parts, ignore_index=True)
    long.to_csv(os.path.join(args.out, "power_combined_long.csv"), index=False, float_format="%.4g")

    # ---- 3. trends
    out_rows = []
    annual_tables = {}
    for r in recs:
        for p in ("T2M", "PRECTOTCORR"):
            yrs, arr = series_array(r["df"], p)
            arr3 = np.where(arr == r["fill"], np.nan, arr)[None, :, :]
            t0 = (yrs[0] + yrs[-1]) / 2
            for sname, months in SEASONS.items():
                y = season_series(arr3, months)[0]
                res = mann_kendall(yrs - t0, y, 8)
                if res is None:
                    continue
                level = float(np.nanmean(y))
                dec_span = (np.isfinite(y).sum() - 1) / 10
                row = dict(النقطة=r["ar"], lat=r["lat"], lon=r["lon"], parameter=p, المتغير=PNAME[p],
                           الموسم=SEASON_AR[sname], season=sname, المتوسط=level,
                           الاتجاه_لكل_عقد=res["slope"] * 10, p_خام=res["p_raw"], p_معدّل=res["p"],
                           ارتباط_ذاتي_lag1=res["r1"], سنين=res["n"],
                           التغير_التقريبي_خلال_الفترة=res["slope"] * 10 * dec_span)
                out_rows.append(row)
                if sname == "annual":
                    annual_tables.setdefault(p, {})[r["ar"]] = pd.Series(y, index=yrs)
    tr = pd.DataFrame(out_rows)
    tr["q_FDR"] = np.nan
    for p in tr["parameter"].unique():
        m = tr["parameter"] == p
        tr.loc[m, "q_FDR"] = bh_qvalues(tr.loc[m, "p_معدّل"].to_numpy())
    tr["معنوي"] = np.where(tr["q_FDR"] < args.alpha, "نعم", "لا")
    tr["الوحدة"] = tr["parameter"].map(lambda x: UNITS[x] + " لكل عقد")
    tr.to_csv(os.path.join(args.out, "power_trends.csv"), index=False, float_format="%.5g")

    show = tr[tr["season"] == "annual"][["النقطة", "parameter", "المتوسط", "الاتجاه_لكل_عقد", "p_معدّل", "q_FDR", "معنوي"]]
    print("\nAnnual trends:\n" + show.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\nSeasonal T2M trends (significant after FDR):")
    s = tr[(tr["parameter"] == "T2M") & (tr["season"] != "annual") & (tr["معنوي"] == "نعم")]
    print(s[["النقطة", "season", "الاتجاه_لكل_عقد", "q_FDR"]].to_string(index=False, float_format=lambda x: f"{x:.4g}")
          if len(s) else "none")

    # ---- Excel
    notes = pd.DataFrame({"ملاحظات": [
        "المصدر: NASA POWER (MERRA-2)، شهري 2001–2025، بيانات حقيقية نزّلها الفريق من POWER Data Access Viewer.",
        "MERRA-2 نموذج إعادة تحليل (reanalysis) مش محطات رصد، ودقته حوالي 0.5° × 0.625° (~50 كم). النقطة = متوسط خلية.",
        "الترند بـ Mann-Kendall + Sen's slope على سلاسل سنوية وموسمية بعد إزالة الموسمية. p معدّل للارتباط الذاتي.",
        "q = تصحيح Benjamini-Hochberg عبر كل اختبارات نفس المتغير (5 نقط × 5 مواسم) علشان ما نطلعش 'اكتشافات' بالصدفة.",
        "الاتجاه المعنوي بيثبت إن فيه تغيير، مش سببه. 25 سنة فترة قصيرة، وسنة شاذة زي 2020 (أمطار) بتأثر على الأمطار.",
        "مفيش NDVI ولا عمران هنا: الملفات دي حرارة وأمطار بس.",
        "الفحص الكامل للملفات في شيت 'فحص الملفات'.",
    ]})
    with pd.ExcelWriter(os.path.join(args.out, "power_real_results.xlsx"), engine="openpyxl") as w:
        notes.to_excel(w, "اقرأني", index=False)
        validation.to_excel(w, "فحص الملفات", index=False)
        tr[["النقطة", "المتغير", "الموسم", "المتوسط", "الاتجاه_لكل_عقد", "الوحدة", "التغير_التقريبي_خلال_الفترة",
            "p_خام", "p_معدّل", "q_FDR", "معنوي", "ارتباط_ذاتي_lag1", "سنين"]].round(5).to_excel(w, "الاتجاهات", index=False)
        for p, d in annual_tables.items():
            pd.DataFrame(d).round(3).rename_axis("السنة").to_excel(w, f"سنوي {p}")
        for ws in w.book.worksheets:
            ws.sheet_view.rightToLeft = True
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(70, max(10, width + 2))
    print(f"\nWrote {args.out}/power_combined_long.csv, power_trends.csv, power_real_results.xlsx")


if __name__ == "__main__":
    main()

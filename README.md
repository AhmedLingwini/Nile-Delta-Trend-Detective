# Nile Delta Trend Detective

NASA Space Apps Challenge: *Be An Earth System Trend Detective*.
Is the Nile Delta getting greener, and can we trust the trend? We analyse 25 years (2001–2025) of NDVI for 468 grid cells (0.1° ≈ 10 km), together with climate (NASA POWER) and night lights (VIIRS), using Theil–Sen slopes, Mann–Kendall tests, autocorrelation correction and FDR.

**مشروع لتحدي NASA Space Apps.** هل دلتا النيل بتخضر؟ وهل نصدق الترند؟ حللنا 25 سنة من NDVI لـ 468 خلية، مع المناخ والإضاءة الليلية، مع اختبارات معنوية وتصحيحات علشان نعرف مين نصدقه.

## Main findings / أهم النتايج

- Regional mean NDVI rose from 0.309 (2001) to 0.339 (2025); Mann-Kendall p ≈ 5.6e-07.
- 382 of 468 cells keep a significant trend after autocorrelation correction and FDR: 204 greening, 178 browning (of these, high confidence: 202 greening, 165 browning).
- Near four of five cities (Alexandria, Mansoura, Damietta, Cairo) NDVI declines within 15 km; near Ismailia it rises.
- Temperature rises at all five POWER points (+0.29 to +0.43 °C per decade, significant after FDR). Rainfall is significant after FDR only at Ismailia.
- Night lights do not rise in any city (a significant decline only in Cairo), so they do not confirm an urban-sprawl story.
- Correlations with climate do not prove causes (NDVI vs rainfall: raw r = +0.43, de-trended r = −0.06).

## Data sources / المصادر

| Data | Source |
|---|---|
| NDVI | NASA Terra MODIS monthly vegetation indices (MOD13C2 / MOD13A3), 2001–2025, aggregated to a 0.1° grid |
| Temperature (T2M), precipitation (PRECTOTCORR) | NASA POWER monthly point data (MERRA-2 reanalysis, about 50 km), 5 points, 2001–2025 |
| Night lights | Suomi NPP VIIRS Day/Night Band monthly composites (NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG, band avg_rad, nW/cm²/sr), downloaded from Google Earth Engine with `earth_engine_nightlights.js`: mean within 15 km of each city, April 2012–December 2025 |

POWER values come from a reanalysis model, not from weather stations.

## Repository layout

```
trend_detective.py          NDVI trends per cell, confidence levels, regional trend, figures
combine_power.py            validates and analyses the POWER files (needs trend_analysis.py, config.py)
trend_analysis.py           Theil-Sen + Mann-Kendall, BH-FDR, seasonal means
config.py                   season definitions
city_ndvi_trends.py         NDVI trend within 15 km of each city
city_nightlights_trends.py  night-light trends per city
earth_engine_nightlights.js Earth Engine script that exports the night-light CSVs
requirements.txt
data/ndvi/                  ndvi_2001_2005.csv ... ndvi_2021_2025.csv   (cell_id, lon, lat, date, mean)
data/power/                 POWER_Point_Monthly_*.csv
data/nightlights/           nightlights_<City>.csv                       (date, mean)
results/                    outputs
```

## How to run

```
pip install -r requirements.txt

python trend_detective.py --ndvi_dir data/ndvi --power_dir data/power --out results
python combine_power.py --data data/power --out results/power
python city_ndvi_trends.py data/ndvi
python city_nightlights_trends.py data/nightlights
```

## Methods / الطريقة

1. Seasonality is removed by subtracting each calendar month's mean.
2. Trend size: Theil–Sen slope. Significance: Mann–Kendall.
3. Autocorrelation: monthly values resemble the previous month, so p-values are corrected. The three scripts use slightly different corrections:
   - `trend_detective.py`: Hamed–Rao correction on monthly anomalies.
   - `combine_power.py`: variance of S inflated by (1 + r1) / (1 − r1) when lag-1 autocorrelation r1 > 0 (annual and seasonal series).
   - `city_*.py`: rank-based Hamed–Rao correction (never below 1) on monthly anomalies.
4. Multiple testing: Benjamini–Hochberg FDR across all cells (or across all POWER tests of one variable).
5. Confidence levels per cell: high = significant after the corrections (FDR q < 0.05) and the annual-means test agrees (same direction, p < 0.05); medium = significant after the corrections only; low = otherwise.
6. Visual check: before/after imagery in Google Earth for 3 cells.

## Limitations / الحدود

- A cell is about 10 km and mixes farmland, buildings and water.
- Climate comes from 5 points of a reanalysis model, so it cannot explain each cell.
- 25 years is short, and sensor changes can affect part of the trend.
- Significant trends show that something changed, not why.
- Night-light trends depend on the product and the radius used (here a 15 km circle around each city point).
- Night-light brightness is not built-up area, and the DNB sensor (about 0.5–0.9 µm) may miss blue light, so LED conversion can look like a decline.
- Linear trends only.

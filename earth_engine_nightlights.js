// Earth Engine (code.earthengine.google.com) script that exports the night-light CSVs used in this project.
// Product: NOAA VIIRS DNB monthly composites (Suomi NPP), band avg_rad (nW/cm^2/sr).
// For each city: mean radiance inside a 15 km circle around the city point, pixels with cf_cvg == 0 masked out.
// Run, then start the 5 export tasks in the Tasks tab. Each CSV has columns: date (YYYY-MM), mean.
// The collection starts in 2012-04, so the first month in the files is 2012-04.

var COLLECTION = 'NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG';   // alternative: .../VCMSLCFG (stray-light corrected, starts 2014-01)
var BAND = 'avg_rad';
var START = '2012-01-01';
var END = '2026-01-01';          // exclusive: last month = December 2025
var RADIUS_M = 15000;            // buffer radius in metres

var CITIES = {                   // [lon, lat]
  Alexandria: [29.92, 31.20],
  Mansoura:   [31.38, 31.04],
  Damietta:   [31.82, 31.42],
  Cairo:      [31.24, 30.04],
  Ismailia:   [32.27, 30.59]
};

var col = ee.ImageCollection(COLLECTION).filterDate(START, END);
print('images', col.size());
print('first month', ee.Date(col.aggregate_min('system:time_start')));
print('last month', ee.Date(col.aggregate_max('system:time_start')));

Object.keys(CITIES).forEach(function (name) {
  var geom = ee.Geometry.Point(CITIES[name]).buffer(RADIUS_M);
  var fc = col.map(function (img) {
    var rad = img.select(BAND).updateMask(img.select('cf_cvg').gt(0));
    var v = rad.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: geom, scale: 500, maxPixels: 1e9
    }).get(BAND);
    return ee.Feature(null, {date: img.date().format('YYYY-MM'), mean: v});
  }).filter(ee.Filter.notNull(['mean']));
  Export.table.toDrive({
    collection: fc, description: 'nightlights_' + name,
    fileFormat: 'CSV', selectors: ['date', 'mean']
  });
});

// Google Earth Engine script: reference ET products over the Karun-Karkheh-Marzi-e-Gharb basins.
//
// 1. Basin-mean monthly series (mm/d): one value per monthly window. Windows start at the first image
//    of each collection and advance by one month, as in data/products/et_products_basin_monthly.csv.
//    TerraClimate is exported as monthly total / 30 and MODIS composites as composite total / 8;
//    etcompare/products.py assigns windows to calendar months and rescales TerraClimate.
// 2. Long-term mean ET maps (mm/d) on a 0.01-degree grid aligned with the 0.5-degree CWatM grid
//    (data/products/maps/ET_mean_<name>.tif).
//
// Set BASINS to the uploaded basin polygon (data/basin/3basins.shp) before running.

var BASINS = ee.FeatureCollection('users/YOUR_USER/3basins');
var region = BASINS.geometry();
var TRANSFORM = [0.01, 0, 44.5, 0, -0.01, 37.0];
var EXPORT_REGION = ee.Geometry.Rectangle([44.5, 30.0, 52.5, 37.0], null, false);

function asRate(rate, img) {
  return rate.toFloat().rename('ET').set('system:time_start', img.get('system:time_start'));
}

// number of days in an 8-day MODIS composite (the last composite of a year is shorter)
function compositeDays(img) {
  var d = ee.Date(img.get('system:time_start'));
  var y = d.get('year');
  var daysInYear = ee.Date.fromYMD(ee.Number(y).add(1), 1, 1).difference(ee.Date.fromYMD(y, 1, 1), 'day');
  return ee.Number(8).min(daysInYear.subtract(d.getRelative('day', 'year')));
}

var PRODUCTS = [
  { name: 'GLDAS', id: 'NASA/GLDAS/V021/NOAH/G025/T3H', scale: 27830, start: '2000-01-01', end: '2020-01-01',
    series: function (img) { return asRate(img.select('Evap_tavg').multiply(86400), img); },
    map: function (img) { return asRate(img.select('Evap_tavg').multiply(86400), img); } },
  { name: 'TerraClimate', id: 'IDAHO_EPSCOR/TERRACLIMATE', scale: 4638, start: '1958-01-01', end: '2020-01-01',
    series: function (img) { return asRate(img.select('aet').multiply(0.1).divide(30), img); },
    map: function (img) {
      var d = ee.Date(img.get('system:time_start'));
      return asRate(img.select('aet').multiply(0.1).divide(d.advance(1, 'month').difference(d, 'day')), img);
    } },
  { name: 'MODIS_1km', id: 'MODIS/NTSG/MOD16A2/105', scale: 1000, start: '2000-01-01', end: '2015-01-01',
    series: function (img) { var e = img.select('ET'); return asRate(e.updateMask(e.lt(32700)).multiply(0.1).divide(8), img); },
    map: function (img) { var e = img.select('ET'); return asRate(e.updateMask(e.lt(32700)).multiply(0.1).divide(compositeDays(img)), img); } },
  { name: 'MODIS_500m', id: 'MODIS/006/MOD16A2', scale: 500, start: '2001-01-01', end: '2020-01-01',
    series: function (img) { var e = img.select('ET'); return asRate(e.updateMask(e.lt(32700)).multiply(0.1).divide(8), img); },
    map: function (img) { var e = img.select('ET'); return asRate(e.updateMask(e.lt(32700)).multiply(0.1).divide(compositeDays(img)), img); } },
  { name: 'PML_V2', id: 'CAS/IGSNRR/PML/V2_v017', scale: 500, start: '2000-02-26', end: '2020-01-01',
    series: function (img) { return asRate(img.select(['Ec', 'Es', 'Ei']).reduce(ee.Reducer.sum()), img); },
    map: function (img) { return asRate(img.select(['Ec', 'Es', 'Ei']).reduce(ee.Reducer.sum()), img); } },
  { name: 'SMAP_L4', id: 'NASA/SMAP/SPL4SMGP/007', scale: 11000, start: '2015-03-31', end: '2020-01-01',
    series: function (img) { return asRate(img.select('land_evapotranspiration_flux').multiply(86400), img); },
    map: function (img) { return asRate(img.select('land_evapotranspiration_flux').multiply(86400), img); } }
];
var MAP_NAME = { GLDAS: 'GLDAS', TerraClimate: 'TerraClimate', MODIS_1km: 'MODIS_1km', MODIS_500m: 'MODIS_500m',
                 PML_V2: 'PML', SMAP_L4: 'SMAP' };
var MAP_START = { TerraClimate: '2000-01-01', PML_V2: '2000-03-01', SMAP_L4: '2015-04-01' };

// ---------- 1. basin-mean monthly series ----------
var series = PRODUCTS.map(function (p) {
  var col = ee.ImageCollection(p.id).filterDate(p.start, p.end).map(p.series);
  var first = ee.Date(col.first().get('system:time_start'));
  var n = ee.Date(p.end).difference(first, 'month').floor();
  return ee.FeatureCollection(ee.List.sequence(0, n.subtract(1)).map(function (k) {
    var s = first.advance(k, 'month');
    var c = col.filterDate(s, s.advance(1, 'month'));
    var v = ee.Algorithms.If(c.size().gt(0),
      c.mean().reduceRegion({ reducer: ee.Reducer.mean(), geometry: region, scale: p.scale, maxPixels: 1e13, tileScale: 4 }).get('ET'),
      null);
    return ee.Feature(null, { dataset: p.name, date: s.format('YYYY-MM-dd'), ET_mm_d: v, n_images: c.size() });
  }));
}).reduce(function (a, b) { return a.merge(b); });

Export.table.toDrive({
  collection: series, description: 'et_products_basin_monthly', folder: 'GEE_ET',
  selectors: ['dataset', 'date', 'ET_mm_d', 'n_images'], fileFormat: 'CSV'
});

// ---------- 2. long-term mean maps over the overlap with the CWatM retrospective run ----------
PRODUCTS.forEach(function (p) {
  var start = MAP_START[p.name] || p.start;
  var mean = ee.ImageCollection(p.id).filterDate(start, p.end).map(p.map).mean();
  Export.image.toDrive({
    image: mean.toFloat(), description: 'ET_mean_' + MAP_NAME[p.name], folder: 'GEE_ET',
    fileNamePrefix: 'ET_mean_' + MAP_NAME[p.name], region: EXPORT_REGION, crs: 'EPSG:4326',
    crsTransform: TRANSFORM, maxPixels: 1e10, fileFormat: 'GeoTIFF'
  });
});

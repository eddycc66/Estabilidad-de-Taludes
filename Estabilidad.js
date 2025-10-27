// =============================================================================
// ANÁLISIS MULTICRITERIO DE ESTABILIDAD DE TALUDES - MUNICIPIO LICOMA
// Evaluación de riesgos geotécnicos para planificación de infraestructura
// =============================================================================
// Docente: MSc. Edwin Calle Condori

// Cargar el área de estudio desde Assets
var licoma = ee.FeatureCollection('projects/eddycc66/assets/licoma');
var areaEstudio = licoma.geometry();
var fechaInicio = '2020-01-01';
var fechaFin = '2024-01-01';

Map.centerObject(areaEstudio, 11);
Map.addLayer(areaEstudio, {color: '0000FF'}, 'Área de Estudio - Licoma');

// Calcular estadísticas básicas del área
var areaHectareas = areaEstudio.area().divide(10000);
print('Área total de la cuenca:', areaHectareas.round(), 'hectáreas');

// 1. CARGA Y PROCESAMIENTO DE DATOS MULTIFUENTE
// Modelo Digital de Elevación (SRTM - 30m)
var dem = ee.Image('USGS/SRTMGL1_003').select('elevation').clip(areaEstudio);

// Calcular parámetros morfométricos
var pendiente = ee.Terrain.slope(dem);
var aspecto = ee.Terrain.aspect(dem);

// Calcular curvatura del terreno
function calcularCurvatura(elevacion) {
  var kernelX = ee.Kernel.fixed(3, 3, [
    [-1, 0, 1],
    [-2, 0, 2],
    [-1, 0, 1]
  ]);
  
  var kernelY = ee.Kernel.fixed(3, 3, [
    [-1, -2, -1],
    [0, 0, 0],
    [1, 2, 1]
  ]);
  
  var dzdx = elevacion.convolve(kernelX);
  var dzdy = elevacion.convolve(kernelY);
  var d2zdx2 = dzdx.convolve(kernelX);
  var d2zdy2 = dzdy.convolve(kernelY);
  var curvatura = d2zdx2.add(d2zdy2).abs();
  
  return curvatura.rename('curvatura');
}

var curvatura = calcularCurvatura(dem);

// 2. ANÁLISIS DE DEFORMACIÓN CON SENTINEL-1
print('Procesando datos Sentinel-1 para análisis de deformación...');

var coleccionSentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(areaEstudio)
  .filterDate(fechaInicio, fechaFin)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV');

var desviacionBackscatter = coleccionSentinel1.reduce(ee.Reducer.stdDev());
var cambiosSuperficiales = desviacionBackscatter.select('VV_stdDev');
var cambiosNorm = cambiosSuperficiales.divide(0.3).min(1).max(0);

// 3. ANÁLISIS DE HUMEDAD Y VEGETACIÓN
var coleccionSentinel2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(areaEstudio)
  .filterDate('2023-01-01', '2023-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var imagenCompuesta = coleccionSentinel2.median().clip(areaEstudio);

var ndwi = imagenCompuesta.normalizedDifference(['B3', 'B8']).rename('NDWI');
var ndvi = imagenCompuesta.normalizedDifference(['B8', 'B4']).rename('NDVI');

// 4. ANÁLISIS HIDROLÓGICO - ÍNDICE DE POSICIÓN TOPOGRÁFICA (TPI)
var tpi = dem.subtract(dem.focal_mean(100, 'circle', 'meters'));
var tpiNorm = tpi.abs().divide(50).min(1).max(0);

// 5. MODELO DE SUSCEPTIBILIDAD MULTICRITERIO AVANZADO
var pendienteNorm = pendiente.divide(50).min(1).max(0);
var curvaturaNorm = curvatura.divide(0.2).min(1).max(0);
var humedadNorm = ndwi.add(0.3).divide(1.3).min(1).max(0);
var tpiNormAjustado = tpi.abs().divide(40).min(1).max(0);

// Pesos según metodología AHP
var pesoPendiente = 0.35;
var pesoDeformacion = 0.25;
var pesoCurvatura = 0.20;
var pesoHumedad = 0.10;
var pesoTPI = 0.10;

// Cálculo del Índice de Susceptibilidad Compuesto
var indiceSusceptibilidad = pendienteNorm.multiply(pesoPendiente)
  .add(cambiosNorm.multiply(pesoDeformacion))
  .add(curvaturaNorm.multiply(pesoCurvatura))
  .add(humedadNorm.multiply(pesoHumedad))
  .add(tpiNormAjustado.multiply(pesoTPI));

// 6. CLASIFICACIÓN DE NIVELES DE RIESGO
var clasesRiesgo = ee.Image(1)
  .where(indiceSusceptibilidad.gte(0.2).and(indiceSusceptibilidad.lt(0.4)), 2)
  .where(indiceSusceptibilidad.gte(0.4).and(indiceSusceptibilidad.lt(0.6)), 3)
  .where(indiceSusceptibilidad.gte(0.6).and(indiceSusceptibilidad.lt(0.8)), 4)
  .where(indiceSusceptibilidad.gte(0.8), 5)
  .rename('riesgo')
  .clip(areaEstudio);

// 7. ANÁLISIS ESTADÍSTICO COMPLETO
print('=== ANÁLISIS GEOESTADÍSTICO - CUENCA LICOMA ===');

// Estadísticas del DEM
var estadisticasDEM = dem.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true).combine(ee.Reducer.stdDev(), '', true),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
});

print('=== ESTADÍSTICAS DE ELEVACIÓN ===');
print('Elevación mínima:', estadisticasDEM.get('elevation_min'), 'msnm');
print('Elevación máxima:', estadisticasDEM.get('elevation_max'), 'msnm');
print('Elevación promedio:', estadisticasDEM.get('elevation_mean'), 'msnm');
print('Desviación estándar:', estadisticasDEM.get('elevation_stdDev'), 'msnm');

// Estadísticas de pendiente
var estadisticasPendiente = pendiente.reduceRegion({
  reducer: ee.Reducer.percentile([5, 25, 50, 75, 95], ['p5', 'p25', 'p50', 'p75', 'p95'])
    .combine(ee.Reducer.mean(), '', true),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
});

print('=== ESTADÍSTICAS DE PENDIENTE ===');
print('Pendiente promedio:', estadisticasPendiente.get('slope_mean'), '°');
print('Percentil 5:', estadisticasPendiente.get('slope_p5'), '°');
print('Percentil 25:', estadisticasPendiente.get('slope_p25'), '°');
print('Mediana (P50):', estadisticasPendiente.get('slope_p50'), '°');
print('Percentil 75:', estadisticasPendiente.get('slope_p75'), '°');
print('Percentil 95:', estadisticasPendiente.get('slope_p95'), '°');

// 8. VISUALIZACIÓN PROFESIONAL
var paletaRiesgo = ['#1a9641', '#a6d96a', '#ffffbf', '#fdae61', '#d7191c'];
var nombresClases = ['Muy Bajo', 'Bajo', 'Medio', 'Alto', 'Muy Alto'];

Map.addLayer(clasesRiesgo, 
  {min: 1, max: 5, palette: paletaRiesgo, opacity: 0.8}, 
  'Índice de Susceptibilidad a Deslizamientos');

Map.addLayer(dem, {min: 500, max: 4000, palette: ['blue', 'green', 'brown', 'white']}, 'Modelo Digital de Elevación', false);
Map.addLayer(pendiente, {min: 0, max: 50, palette: ['green', 'yellow', 'red']}, 'Pendiente (°)', false);
Map.addLayer(curvatura, {min: 0, max: 0.2, palette: ['white', 'red']}, 'Curvatura del Terreno', false);
Map.addLayer(ndwi, {min: -0.3, max: 0.5, palette: ['brown', 'yellow', 'blue']}, 'Índice de Humedad (NDWI)', false);

// 9. GENERACIÓN DE REPORTES ESTADÍSTICOS DETALLADOS
print('=== RESUMEN EJECUTIVO - ANÁLISIS DE RIESGOS ===');

// Calcular distribución de áreas por clase de riesgo
var pixelArea = ee.Image.pixelArea().divide(10000); // m² a hectáreas

// Calcular área por clase individualmente
var area1 = pixelArea.updateMask(clasesRiesgo.eq(1)).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
}).get('area');

var area2 = pixelArea.updateMask(clasesRiesgo.eq(2)).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
}).get('area');

var area3 = pixelArea.updateMask(clasesRiesgo.eq(3)).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
}).get('area');

var area4 = pixelArea.updateMask(clasesRiesgo.eq(4)).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
}).get('area');

var area5 = pixelArea.updateMask(clasesRiesgo.eq(5)).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
}).get('area');

// Manejar valores nulos
area1 = ee.Number(ee.Algorithms.If(area1, area1, 0));
area2 = ee.Number(ee.Algorithms.If(area2, area2, 0));
area3 = ee.Number(ee.Algorithms.If(area3, area3, 0));
area4 = ee.Number(ee.Algorithms.If(area4, area4, 0));
area5 = ee.Number(ee.Algorithms.If(area5, area5, 0));

print('• Muy Bajo:', area1.round(), 'ha');
print('• Bajo:', area2.round(), 'ha');
print('• Medio:', area3.round(), 'ha');
print('• Alto:', area4.round(), 'ha');
print('• Muy Alto:', area5.round(), 'ha');

// Calcular área de alto riesgo
var areaAltoRiesgo = area4.add(area5);
var porcentajeAltoRiesgo = areaAltoRiesgo.divide(areaHectareas).multiply(100);

print('=== INDICADORES CRÍTICOS ===');
print('Área total de alto y muy alto riesgo:', areaAltoRiesgo.round(), 'ha');
print('Porcentaje de área crítica:', porcentajeAltoRiesgo.round(), '%');

// Análisis por rangos de pendiente
var areaPendienteAlta = pixelArea.updateMask(pendiente.gt(30)).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: areaEstudio,
  scale: 30,
  maxPixels: 1e10
}).get('area');

areaPendienteAlta = ee.Number(ee.Algorithms.If(areaPendienteAlta, areaPendienteAlta, 0));
var porcentajePendienteAlta = areaPendienteAlta.divide(areaHectareas).multiply(100);

print('Área con pendiente >30°:', areaPendienteAlta.round(), 'ha');
print('Porcentaje de área con pendiente >30°:', porcentajePendienteAlta.round(), '%');

// 10. GENERACIÓN DE GRÁFICOS
var chartHistograma = ui.Chart.image.histogram({
  image: indiceSusceptibilidad,
  region: areaEstudio,
  scale: 30,
  maxPixels: 1e9,
  minBucketWidth: 0.05
}).setOptions({
  title: 'Distribución del Índice de Susceptibilidad - Cuenca Licoma',
  hAxis: {title: 'Valor del Índice de Susceptibilidad'},
  vAxis: {title: 'Frecuencia (número de píxeles)'},
  colors: ['#1a9641'],
  legend: {position: 'none'}
});

print(chartHistograma);

// Gráfico de distribución de pendientes
var chartPendiente = ui.Chart.image.histogram({
  image: pendiente,
  region: areaEstudio,
  scale: 30,
  maxPixels: 1e9,
  minBucketWidth: 2
}).setOptions({
  title: 'Distribución de Pendientes - Cuenca Licoma',
  hAxis: {title: 'Pendiente (grados)'},
  vAxis: {title: 'Frecuencia (número de píxeles)'},
  colors: ['#d95f0e'],
  legend: {position: 'none'}
});

print(chartPendiente);

// 11. EXPORTACIÓN DE RESULTADOS
Export.image.toDrive({
  image: clasesRiesgo.toByte(),
  description: 'Mapa_Susceptibilidad_Licoma',
  scale: 30,
  region: areaEstudio,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  folder: 'GEE_Exports'
});

Export.image.toDrive({
  image: indiceSusceptibilidad.float(),
  description: 'Indice_Susceptibilidad_Continuo_Licoma',
  scale: 30,
  region: areaEstudio,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  folder: 'GEE_Exports'
});

Export.image.toDrive({
  image: dem,
  description: 'DEM_Licoma',
  scale: 30,
  region: areaEstudio,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  folder: 'GEE_Exports'
});

Export.image.toDrive({
  image: pendiente.float(),
  description: 'Pendiente_Licoma',
  scale: 30,
  region: areaEstudio,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  folder: 'GEE_Exports'
});

Export.image.toDrive({
  image: curvatura.float(),
  description: 'Curvatura_Licoma',
  scale: 30,
  region: areaEstudio,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  folder: 'GEE_Exports'
});

// Exportar estadísticas
var estadisticasExport = ee.FeatureCollection([
  ee.Feature(null, {
    'Area_Total_ha': areaHectareas,
    'Elevacion_Minima': estadisticasDEM.get('elevation_min'),
    'Elevacion_Maxima': estadisticasDEM.get('elevation_max'),
    'Elevacion_Promedio': estadisticasDEM.get('elevation_mean'),
    'Pendiente_Promedio': estadisticasPendiente.get('slope_mean'),
    'Pendiente_P50': estadisticasPendiente.get('slope_p50'),
    'Pendiente_P95': estadisticasPendiente.get('slope_p95'),
    'Area_Muy_Bajo_Riesgo': area1,
    'Area_Bajo_Riesgo': area2,
    'Area_Medio_Riesgo': area3,
    'Area_Alto_Riesgo': area4,
    'Area_Muy_Alto_Riesgo': area5,
    'Area_Alto_Riesgo_Total': areaAltoRiesgo,
    'Porcentaje_Area_Critica': porcentajeAltoRiesgo,
    'Area_Pendiente_Mayor_30': areaPendienteAlta,
    'Porcentaje_Pendiente_Mayor_30': porcentajePendienteAlta
  })
]);

Export.table.toDrive({
  collection: estadisticasExport,
  description: 'Estadisticas_Detalladas_Licoma',
  fileFormat: 'CSV',
  folder: 'GEE_Exports'
});

// 12. LEYENDA INTERACTIVA
function crearLeyenda() {
  var leyenda = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '8px 15px',
      backgroundColor: 'white',
      border: '1px solid #ccc',
      borderRadius: '5px'
    }
  });

  var tituloLeyenda = ui.Label({
    value: 'NIVELES DE RIESGO GEOTÉCNICO',
    style: {
      fontWeight: 'bold',
      fontSize: '14px',
      margin: '0 0 8px 0',
      textAlign: 'center'
    }
  });
  
  leyenda.add(tituloLeyenda);

  for (var i = 4; i >= 0; i--) {
    var colorBox = ui.Label({
      style: {
        backgroundColor: paletaRiesgo[i],
        padding: '10px',
        margin: '2px',
        width: '20px',
        height: '20px',
        border: '1px solid #333'
      }
    });
    
    var descripcion = ui.Label({
      value: nombresClases[i],
      style: {margin: '0 0 0 8px', fontSize: '12px'}
    });
    
    var filaLeyenda = ui.Panel({
      widgets: [colorBox, descripcion],
      layout: ui.Panel.Layout.Flow('horizontal')
    });
    
    leyenda.add(filaLeyenda);
  }
  
  return leyenda;
}

Map.add(crearLeyenda());

// 13. TABLA RESUMEN FINAL
print('=== TABLA RESUMEN - PARÁMETROS DEL MODELO ===');
print('• Peso Pendiente: ' + (pesoPendiente * 100) + '%');
print('• Peso Deformación: ' + (pesoDeformacion * 100) + '%');
print('• Peso Curvatura: ' + (pesoCurvatura * 100) + '%');
print('• Peso Humedad: ' + (pesoHumedad * 100) + '%');
print('• Peso TPI: ' + (pesoTPI * 100) + '%');

print('=== RECOMENDACIONES DE INGENIERÍA ===');
print('• Zonas VERDES: Aptas para construcción sin restricciones');
print('• Zonas AMARILLAS: Requieren estudios geotécnicos básicos');
print('• Zonas NARANJAS: Necesitan obras de contención moderadas');
print('• Zonas ROJAS: Evitar construcción o diseñar obras mayores');

print('=== ANÁLISIS COMPLETADO ===');
print('Todos los cálculos y exportaciones han sido ejecutados exitosamente.');
print('Revisa la pestaña "Tasks" (ícono de engranaje ⚙) en el panel derecho para iniciar las descargas.');
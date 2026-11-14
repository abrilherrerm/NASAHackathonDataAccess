# Contexto del proyecto — IMERG México

Documento para quien se sume al proyecto (personas o asistentes de IA): qué hace, por qué
está hecho así, cómo está organizado, qué problemas ya se resolvieron y qué sigue.
El README explica *cómo correrlo*; este archivo explica *el porqué*.

_Última actualización: 2026-11-14 (límites INEGI y series zonales)_

---

## 1. Objetivo

Hackathon 2026. Construir una herramienta para explorar la **precipitación diaria en México
durante 5 años** con datos satelitales de NASA (GPM IMERG):

1. Extraer los datos (equivalente a lo que se haría a mano en **Giovanni**), pero automatizado con `earthaccess` + `Harmony`.
2. Recortarlos a México en la nube (Harmony) y guardarlos en **NetCDF con chunks**.
3. Servir **series de tiempo interactivas** con **FastAPI** y un frontend web (mapa + gráficas).

## 2. Flujo de datos

```
NASA CMR / Earthdata Login
        │  earthaccess (login + búsqueda de la colección)
        ▼
Harmony (HOSS) ── recorte con polígono de México (data/mexico.geojson), 1 job por mes
        │
        ▼
data/raw/YYYY-MM/*.nc4        ← 1 archivo por día (~80 KB c/u, ~2.5 MB por mes)
        │  process.py (xarray + dask)
        ▼
data/processed/imerg_mexico_daily_YYYY.nc   ← NetCDF anual, comprimido, chunks (366, 40, 40)
        │  zonal.py + límites INEGI (make_boundaries.py)
data/processed/zonal_series.nc              ← series diarias de México, 32 estados y 2,478 municipios
        │  app/main.py (FastAPI, lectura perezosa)
        ▼
API /api/...  +  frontend app/static/index.html (Leaflet + Plotly)
```

## 3. Decisiones clave (y por qué)

| Decisión | Motivo |
|---|---|
| **Producto `GPM_3IMERGDF_07`** (diario Final), no `GPM_3IMERGHH_07` (30 min) | El HH en 5 años son ~87,600 archivos; el DF son ~1,826 y **ya es** la precipitación total diaria (acumulado de los HH). Mismo resultado, ~48× menos trabajo. |
| Concept ID fijo `C2723754864-GES_DISC`, variable `precipitation` (mm/día) | Verificado en CMR. La versión V06 usaba `precipitationCal`; en V07 es `precipitation`. |
| **Periodo 2020-10-01 → 2025-09-30** | El Final Run tiene ~3.5 meses de latencia; en nov-2026 el último dato en CMR era 2025-09-30. Ajustar `START_DATE`/`END_DATE` en `config.py` cuando haya más. |
| **Recorte con polígono GeoJSON en Harmony** (no solo bbox) | El servicio HOSS de esta colección soporta `Shapefile: GeoJSON`. Fuera del polígono queda `_FillValue`. `process.py` reaplica la máscara (regionmask, Natural Earth) como respaldo. |
| **Chunks de descarga: 1 job de Harmony por mes** | Jobs chicos (~30 granulos), reanudables. Ritmo observado: **~1.5 min/mes → ~1.5 h** para 5 años. |
| **Chunks en disco `(time=366, lat=40, lon=40)`**, un archivo por año | Una serie de tiempo por punto lee 1 bloque por año → respuestas en < 1 s. |
| Librerías JS **incluidas en el repo** (`app/static/vendor/`) | Evita depender de CDNs (la URL de cdnjs de Plotly daba 404). Solo el mapa base necesita internet. |
| Mapa base **Esri World Gray Canvas** | CARTO empezó a exigir API key ("API KEY REQUIRED"). Esri no la pide. |
| Leaflet con **`zoomSnap: 1`** | Con zoom fraccionario aparecían costuras blancas entre teselas del mapa base. |
| Ráster dibujado en **canvas reproyectado a Web Mercator** | Un `imageOverlay` lineal en lat/lon queda desalineado con el mapa base; se calcula la posición Y de cada fila con la fórmula de Mercator. |
| Promedios espaciales **ponderados por cos(lat)** | Las celdas de 0.1° se hacen más chicas hacia el norte. |
| Límites **INEGI, Marco Geoestadístico Encuesta Intercensal 2025 (ed. 2026)** | Oficiales, con claves CVE_ENT/CVE_MUN. Del ZIP nacional (2.74 GB) solo se baja el integrado (~261 MB) por HTTP Range. |
| Simplificación con `shapely.coverage_simplify` en EPSG:6372 | Mantiene bordes compartidos idénticos (sin huecos entre vecinos). Web: estados 1 km (~925 KB, ~160 KB gzip); municipios 250 m por entidad (mediana ~46 KB gzip, Oaxaca ~220 KB). Cálculo: 50 m. |
| Estados siempre visibles; **municipios solo de la entidad elegida** | Los municipios de todo el país pesan ~9 MB; por entidad se cargan bajo demanda y se guardan en memoria del navegador. |
| Contornos **solo de referencia** (no seleccionables con clic) | Decisión de diseño: el clic en el mapa siempre consulta un punto. Estados/municipios se consultan con el menú "Consultar". |
| Promedio por **fracción de celda** (`zonal.py`) | La resolución no cambia (0.1°); cada celda pesa por la fracción cubierta × cos(lat). 797 de 2,478 municipios son más chicos que una celda (~118 km²); con "centro de celda" muchos quedarían con 0 celdas o idénticos a su vecino. |
| **Series precalculadas** por estado/municipio (`zonal_series.nc`) | 2,511 regiones × 1,826 días ≈ 18 MB; consulta en ~15 ms (vs ~250 ms al vuelo). Se calculan en ~6 s con una matriz dispersa de pesos (scipy). |
| GZip en FastAPI | GeoJSON y JSON grandes viajan comprimidos (~5×). |

## 4. Estructura del repo

```
imerg_mexico/
├── config.py              # rutas, producto, periodo, bbox, tamaños de chunk
├── download.py            # Paso 1: login + jobs de Harmony por mes → data/raw/
├── process.py             # Paso 2: une, máscara México, NetCDF anuales → data/processed/
├── make_boundaries.py     # Límites INEGI → app/static/boundaries/ y data/boundaries/inegi_zonal.gpkg
├── zonal.py               # Pesos por fracción de celda + series por estado/municipio
├── make_demo_data.py      # Datos SINTÉTICOS para probar la API sin cuenta
├── app/
│   ├── main.py            # API FastAPI
│   └── static/
│       ├── index.html     # Frontend (Leaflet + Geoman + Plotly)
│       ├── boundaries/    # estados.geojson, municipios/{CVE_ENT}.geojson, catalogo.json
│       └── vendor/        # leaflet 1.9.4, geoman 2.17.0, plotly 2.35.2
├── data/
│   ├── raw/               # (ignorado en git) salidas de Harmony
│   ├── processed/         # NetCDF anuales + zonal_series.nc
│   ├── boundaries/        # (ignorado en git) ZIP del INEGI + inegi_zonal.gpkg
│   └── mexico.geojson     # polígono usado para el recorte en Harmony (se genera solo)
├── requirements.txt
├── .env.example           # plantilla de credenciales
├── README.md              # cómo correrlo
└── CONTEXT.md             # este archivo
```

## 5. API

| Método | Ruta | Qué hace |
|---|---|---|
| GET | `/api/info` | Rango de fechas disponible, nº de días, bbox, resolución |
| GET | `/api/timeseries/point?lat&lon&freq&start&end` | Serie en la celda más cercana |
| GET | `/api/timeseries/region?west&south&east&north&freq` | Promedio en rectángulo (por defecto todo México) |
| GET | `/api/timeseries/entity/{cvegeo}?freq&start&end` | Serie precalculada: estado (`20`), municipio (`20067`) o México (`00`) + `area_km2`, `stats.volume_hm3` |
| POST | `/api/timeseries/polygon` | Body `{geometry: GeoJSON, freq, start, end, name}` → promedio en el polígono + `area_km2`, `n_cells`, `stats.volume_hm3` |
| GET | `/api/map?start&end&step` | Acumulado por celda para pintar el mapa |

`freq`: `D` diaria, `W` semanal, `MS` mensual, `YS` anual (sumas). Swagger en `/docs`.

Detalles: el dataset se abre una sola vez (`lru_cache`), así que **tras correr `process.py`
hay que reiniciar uvicorn** (`--reload` solo vigila archivos `.py`). Áreas dibujadas también se
promedian por fracción de celda (`zonal.polygon_weights`). Volumen: `mm × km² / 1000 = hm³`.

## 6. Frontend

- Mapa: acumulado del rango elegido (rampa azul secuencial, percentil 98 como máximo), contornos de estados con nombre al pasar el cursor, lectura de coordenadas y mm bajo el cursor.
- Límites INEGI: estados siempre; "Ver municipios de [entidad]" carga solo esa entidad (línea punteada, nombre y CVEGEO al pasar el cursor). Atribución INEGI en el mapa.
- Selección: **punto** (clic), **área** (rectángulo/polígono con Geoman; editable/movible), **estado/municipio** (menú "Consultar"; se resalta en naranja y se encuadra su estado), **México completo**.
- Gráfica Plotly (línea si es diaria, barras si es agregada) + tarjetas: total, promedio, máximo y fecha; para áreas también km² y hm³.

## 7. Entorno y credenciales

- Python 3.12 en **`.venv`** (no conda). Si aparece `No module named 'xarray'`, se está usando el Python de conda: usar `.venv/bin/python ...` o `conda deactivate && source .venv/bin/activate`.
- Credenciales de **NASA Earthdata Login** en `.env` (`EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`), que **no se sube** a git. Alternativa: `~/.netrc`. La cuenta debe tener autorizada la app *NASA GESDISC DATA ARCHIVE*.
- `h5py` es obligatorio (sin él `process.py` no puede leer ningún archivo).

## 8. Estado actual y pendientes

**Hecho**
- [x] Pipeline completo descarga → proceso → API → frontend, probado con datos reales.
- [x] Selección por punto, área dibujada y México completo; límites estatales.

**En curso**
- [ ] Descarga completa de los 60 meses (al 14-nov-2026 iba en ~2022-11). Al terminar: `process.py` y reiniciar uvicorn.

- [x] Límites oficiales INEGI (estados + municipios por entidad), menú de consulta, fracción de celda, series precalculadas.

**Pendiente menor**
- [ ] El recorte de Harmony (`data/mexico.geojson`) y la máscara de `process.py` siguen usando Natural Earth. Podría cambiarse a la unión de estados INEGI en una futura re-descarga (no urge: la diferencia es solo en celdas de costa).

**Ideas**
- Anomalías respecto a climatología, días con lluvia > umbral, exportar CSV de la serie.

## 9. Problemas ya resueltos (para no repetirlos)

| Síntoma | Causa | Solución |
|---|---|---|
| `TypeError: 'str' object is not callable` en `concept_id()` | API de earthaccess cambió | Se usa el concept ID fijo en `config.py` |
| `process.py`: "Abriendo 0 archivos" | Faltaba `h5py` | Agregado a `requirements.txt` |
| Meses incompletos tras cortar la descarga | Se saltaba cualquier carpeta con ≥1 archivo | Ahora exige tantos archivos como días del mes |
| Datos sintéticos mezclados con reales | `process.py` solo sobrescribía años con datos | Borra los NetCDF procesados antes de escribir |
| Front: `Can't find variable: Plotly` | URL de cdnjs con 404 | Plotly local en `vendor/` |
| Mapa base con "API KEY REQUIRED" | CARTO exige clave | Cambio a Esri |
| Mapa gris vacío | Comentario `//` a media línea anuló `.fitBounds()` | Comentario en su propia línea |
| Descarga atorada en un mes (sin avanzar) | Descarga HTTP colgada | `Ctrl+C` y volver a correr `download.py` (reanuda y repite el mes incompleto) |

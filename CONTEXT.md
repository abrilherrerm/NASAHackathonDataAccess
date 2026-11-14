# Contexto del proyecto — IMERG México

Documento para quien se sume al proyecto (personas o asistentes de IA): qué hace, por qué
está hecho así, cómo está organizado, qué problemas ya se resolvieron y qué sigue.
El README explica *cómo correrlo*; este archivo explica *el porqué*.

_Última actualización: 2026-11-14_

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
| Límites de estados de **Natural Earth 1:10m** (provisional) | Rápido y ligero (121 KB) para ubicar al usuario. **Se planea reemplazar por INEGI** (ver §8). |

## 4. Estructura del repo

```
imerg_mexico/
├── config.py              # rutas, producto, periodo, bbox, tamaños de chunk
├── download.py            # Paso 1: login + jobs de Harmony por mes → data/raw/
├── process.py             # Paso 2: une, máscara México, NetCDF anuales → data/processed/
├── make_states.py         # Genera app/static/mexico_states.geojson (Natural Earth)
├── make_demo_data.py      # Datos SINTÉTICOS para probar la API sin cuenta
├── app/
│   ├── main.py            # API FastAPI
│   └── static/
│       ├── index.html     # Frontend (Leaflet + Geoman + Plotly)
│       ├── mexico_states.geojson
│       └── vendor/        # leaflet 1.9.4, geoman 2.17.0, plotly 2.35.2
├── data/
│   ├── raw/               # (ignorado en git) salidas de Harmony
│   ├── processed/         # NetCDF anuales
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
| POST | `/api/timeseries/polygon` | Body `{geometry: GeoJSON, freq, start, end, name}` → promedio en el polígono + `area_km2`, `n_cells`, `stats.volume_hm3` |
| GET | `/api/map?start&end&step` | Acumulado por celda para pintar el mapa |

`freq`: `D` diaria, `W` semanal, `MS` mensual, `YS` anual (sumas). Swagger en `/docs`.

Detalles: el dataset se abre una sola vez (`lru_cache`), así que **tras correr `process.py`
hay que reiniciar uvicorn** (`--reload` solo vigila archivos `.py`). Polígonos más chicos que
una celda usan la celda del centroide. Volumen: `mm × km² / 1000 = hm³`.

## 6. Frontend

- Mapa: acumulado del rango elegido (rampa azul secuencial, percentil 98 como máximo), contornos de estados con nombre al pasar el cursor, lectura de coordenadas y mm bajo el cursor.
- Selección: **punto** (clic), **área** (rectángulo/polígono con Geoman; editable/movible), **México completo**.
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

**Siguiente (acordado)**
- [ ] **Límites oficiales del INEGI (Marco Geoestadístico): estados (32) y municipios (2,475)**
  - Reproyectar de EPSG:6372 (Lambert) a EPSG:4326.
  - Simplificar **preservando topología** (mapshaper/TopoJSON) para que no queden huecos; objetivo ~1–3 MB.
  - Municipios visibles solo con zoom cercano; tooltip con nombre y claves `CVE_ENT`/`CVE_MUN`.
  - Atribución "Fuente: INEGI, Marco Geoestadístico" en el mapa.
- [ ] Seleccionar un estado/municipio como área de consulta.
- [ ] **Estadística zonal ponderada por fracción de celda** (p. ej. `exactextract`): muchos municipios son más chicos que una celda IMERG (~11×11 km), y el método actual (centro de celda dentro del polígono) da 0–1 celdas.

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

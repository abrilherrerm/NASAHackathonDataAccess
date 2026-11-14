# IMERG V07 – Precipitación diaria de México (5 años)

Pipeline: **earthaccess** (login + búsqueda CMR) → **Harmony** (recorte de México en la nube, por chunks mensuales)
→ **xarray** (total diario, máscara del polígono de México, NetCDF anuales con chunks) → **FastAPI** (series de tiempo interactivas).

## Credenciales (NASA Earthdata Login)
Opción A — archivo `.env` (recomendado):
```bash
cp .env.example .env    # y edita usuario/contraseña
```
Opción B — no hagas nada: la primera vez `python download.py` te pide usuario/contraseña
en la terminal y los guarda en `~/.netrc` para siempre.

## Uso
```bash
source .venv/bin/activate
python download.py --start 2025-01-01 --end 2025-01-31   # prueba de 1 mes
python download.py                                     # 5 años: 2020-10-01 → 2025-09-30 (60 jobs)
python make_boundaries.py                              # límites INEGI (una vez, ~5 min la 1ª vez)
python process.py                                      # -> NetCDF anuales + series por estado/municipio
uvicorn app.main:app --reload                          # http://127.0.0.1:8000  ·  /docs
```
Sin cuenta, para probar la API: `python make_demo_data.py` (datos SINTÉTICOS; process.py los sobrescribe).

## Producto
**GPM_3IMERGDF_07** (`C2723754864-GES_DISC`): IMERG Final Run **diario**, 0.1°, global (cubre México),
variable `precipitation` en mm/día (= total diario acumulado de los HH de 30 min).
Disponible 1998-01-01 → 2025-09-30 (latencia ~3.5 meses; ajusta `START_DATE`/`END_DATE` en `config.py` cuando haya más).
Harmony (servicio HOSS) lo recorta con el polígono de México (`data/mexico.geojson`).

## Endpoints
- `GET /api/info`
- `GET /api/timeseries/point?lat=19.43&lon=-99.13&freq=D|W|MS|YS&start=&end=`
- `GET /api/timeseries/region?west=&south=&east=&north=&freq=` (promedio ponderado por cos(lat); por defecto todo México)
- `POST /api/timeseries/polygon` — body `{geometry: GeoJSON, freq, start, end, name}`: promedio ponderado por área dentro del polígono + área (km²) y volumen precipitado (hm³)
- `GET /api/timeseries/entity/{cvegeo}?freq=` — serie precalculada de un estado (`20`), municipio (`20067`) o México (`00`), ponderada por fracción de celda
- `GET /api/map?start=2025-06-01&end=2025-09-30&step=2` (acumulado para el mapa)

## Chunks
- **Descarga:** un job de Harmony por mes (~30 archivos diarios) (`REQUEST_CHUNK` en `config.py`); si se corta, al reanudar omite los meses ya descargados.
- **Disco:** NetCDF por año, comprimido, `chunksizes=(366, 40, 40)` → una serie por punto lee un solo bloque por año.

## Frontend (http://127.0.0.1:8000)
- Mapa Leaflet con el acumulado de precipitación y los **32 estados del INEGI** (referencia, nombre al pasar el cursor).
- **Ver municipios de [entidad]**: carga solo los municipios de esa entidad (ocultos por defecto).
- Consultas: **punto** (clic), **área** dibujada (rectángulo/polígono), **estado/municipio** (menú "Consultar"), **México completo**.
- Áreas, estados y municipios se promedian **por fracción de celda** (`zonal.py`); la resolución IMERG (0.1°) no cambia.
- Librerías JS en `app/static/vendor/`; solo el mapa base (Esri) requiere internet.

## Límites (INEGI)
`make_boundaries.py` baja solo el paquete integrado (~261 MB) del Marco Geoestadístico, Encuesta Intercensal 2025
(ed. 2026), y genera `app/static/boundaries/` (estados ~1 km, municipios por entidad ~250 m, catálogo) y
`data/boundaries/inegi_zonal.gpkg` (~50 m, para el cálculo). Fuente: INEGI.

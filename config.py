"""Configuración central del pipeline IMERG diario -> México."""
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "data" / "raw"              # salidas de Harmony (recortadas)
PROCESSED_DIR = BASE_DIR / "data" / "processed"  # NetCDF diarios chunked
MEXICO_GEOJSON = BASE_DIR / "data" / "mexico.geojson"

# --- Producto: IMERG Final Run diario V07 (mm/día), global 0.1° -------------
# Ya trae la precipitación total diaria (acumulado de los HH de 30 min).
SHORT_NAME = "GPM_3IMERGDF"
VERSION = "07"
CONCEPT_ID = "C2723754864-GES_DISC"   # verificado en CMR (nov-2026)
VARIABLE = "precipitation"

# --- Periodo: 5 años. El Final Run tiene ~3.5 meses de latencia; -------------
# al 14-nov-2026 el último dato disponible es 2025-09-30.
START_DATE = "2020-10-01"
END_DATE = "2025-09-30"

# --- México: bbox (para metadatos/API) y polígono GeoJSON (para Harmony) -----
MEXICO_BBOX = {"west": -118.5, "south": 14.4, "east": -86.6, "north": 32.8}

# --- Chunking -----------------------------------------------------------------
REQUEST_CHUNK = "MS"   # un job de Harmony por mes (~30 granulos)
NC_CHUNKS = {"time": 366, "lat": 40, "lon": 40}   # chunks en disco

# --- Límites INEGI y series zonales (estados/municipios) ----------------------
BOUNDARIES_DIR = BASE_DIR / "data" / "boundaries"
ZONAL_GPKG = BOUNDARIES_DIR / "inegi_zonal.gpkg"          # geometrías para el cálculo
ZONAL_SERIES = PROCESSED_DIR / "zonal_series.nc"          # series diarias precalculadas

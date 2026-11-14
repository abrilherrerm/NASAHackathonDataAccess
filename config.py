"""Configuración central del pipeline IMERG diario -> México."""
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "data" / "raw"              # salidas de Harmony (recortadas)
PROCESSED_DIR = BASE_DIR / "data" / "processed"  # NetCDF diarios chunked
MEXICO_GEOJSON = BASE_DIR / "data" / "mexico.geojson"

# --- Productos IMERG diarios V07 (mm/día), global 0.1° ------------------------
# Final: calibrado con pluviómetros (GPCC), ~3.5 meses de latencia  -> referencia
# Late : preliminar, ~14 h de latencia, sin calibración mensual     -> solo donde no hay Final
RUNS = {
    "final": {"short_name": "GPM_3IMERGDF", "concept_id": "C2723754864-GES_DISC", "raw_dir": RAW_DIR},
    "late":  {"short_name": "GPM_3IMERGDL", "concept_id": "C2723754859-GES_DISC",
              "raw_dir": BASE_DIR / "data" / "raw_late"},
}
VERSION = "07"
VARIABLE = "precipitation"

# Compatibilidad con código previo
SHORT_NAME = RUNS["final"]["short_name"]
CONCEPT_ID = RUNS["final"]["concept_id"]

# --- Periodo -------------------------------------------------------------------
# Final hasta su último día disponible; después, Late hasta el último disponible.
# download.py consulta CMR para saber ambas fechas (END_DATE=None -> lo más reciente).
START_DATE = "2014-06-01"   # inicio de la era GPM (mejores sensores)
END_DATE = None

# --- México: bbox (para metadatos/API) y polígono GeoJSON (para Harmony) -----
MEXICO_BBOX = {"west": -118.5, "south": 14.4, "east": -86.6, "north": 32.8}

# --- Chunking -----------------------------------------------------------------
REQUEST_CHUNK = "MS"   # un job de Harmony por mes (~30 granulos)
NC_CHUNKS = {"time": 366, "lat": 40, "lon": 40}   # chunks en disco

# --- Límites INEGI y series zonales (estados/municipios) ----------------------
BOUNDARIES_DIR = BASE_DIR / "data" / "boundaries"
ZONAL_GPKG = BOUNDARIES_DIR / "inegi_zonal.gpkg"          # geometrías para el cálculo
ZONAL_SERIES = PROCESSED_DIR / "zonal_series.nc"          # series diarias precalculadas

# --- IMERG Final V07 MENSUAL (climatologías, SPI, tendencias de largo plazo) ----
MONTHLY_CONCEPT_ID = "C2723754851-GES_DISC"   # GPM_3IMERGM v07, verificado en CMR (nov-2026): 1998-01 → 2025-09
MONTHLY_START = "1998-01-01"
RAW_MONTHLY_DIR = BASE_DIR / "data" / "raw_monthly"
MONTHLY_FILE = PROCESSED_DIR / "imerg_mexico_monthly.nc"
ZONAL_MONTHLY = PROCESSED_DIR / "zonal_monthly.nc"

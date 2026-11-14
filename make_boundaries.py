"""Genera los límites oficiales del INEGI (Marco Geoestadístico, Encuesta Intercensal 2025)
para el mapa y para el cálculo zonal.

    python make_boundaries.py            # descarga (si hace falta) y genera todo
    python make_boundaries.py ruta.zip   # usa un MG_Integrado_*.zip ya descargado

Del ZIP nacional (2.74 GB) solo se baja el paquete integrado (~261 MB) usando HTTP Range,
y de él solo se extraen 00ent (estados) y 00mun (municipios).

Salidas:
  app/static/boundaries/estados.geojson          32 estados, simplificado a ~1 km (web)
  app/static/boundaries/municipios/{CVE_ENT}.geojson   municipios por entidad, ~250 m (web)
  app/static/boundaries/catalogo.json            estados y municipios para los menús
  data/boundaries/inegi_zonal.gpkg               geometrías ~50 m para el cálculo zonal
"""
import io
import json
import struct
import sys
import urllib.request
import zipfile
from pathlib import Path

import geopandas as gpd
import pandas as pd
import shapely

import config

MG_URL = ("https://www.inegi.org.mx/contenidos/productos/prod_serv/contenidos/espanol/bvinegi/"
          "productos/geografia/marcogeo/794551196649_s.zip")
INNER = "MG_Integrado_Encuesta_Intercensal_2025.zip"
SOURCE = "INEGI, Marco Geoestadístico, Encuesta Intercensal 2025 (ed. 2026)"

CACHE_ZIP = config.BOUNDARIES_DIR / INNER
WEB_DIR = config.BASE_DIR / "app" / "static" / "boundaries"
TOL_ESTADOS_M, TOL_MUNICIPIOS_M, TOL_ZONAL_M = 1000, 250, 50


# ---------------------------------------------------------------- descarga parcial
class _HTTPRangeFile(io.RawIOBase):
    """Archivo remoto de solo lectura que pide únicamente los bytes que se leen."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        r = urllib.request.urlopen(self._req(0, 0))
        self.size = int(r.headers["Content-Range"].split("/")[-1])

    def _req(self, a, b):
        return urllib.request.Request(self.url, headers={"User-Agent": "Mozilla/5.0", "Range": f"bytes={a}-{b}"})

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        data = urllib.request.urlopen(self._req(self.pos, self.pos + n - 1), timeout=120).read()
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def download_integrated() -> Path:
    if CACHE_ZIP.exists():
        return CACHE_ZIP
    print(f"Descargando solo {INNER} (~261 MB) del Marco Geoestadístico…")
    CACHE_ZIP.parent.mkdir(parents=True, exist_ok=True)
    outer = zipfile.ZipFile(io.BufferedReader(_HTTPRangeFile(MG_URL), buffer_size=1 << 20))
    tmp = CACHE_ZIP.with_suffix(".part")
    with outer.open(INNER) as src, open(tmp, "wb") as dst:
        while chunk := src.read(1 << 22):
            dst.write(chunk)
    tmp.rename(CACHE_ZIP)
    return CACHE_ZIP


def read_layer(zip_path: Path, name: str) -> gpd.GeoDataFrame:
    return gpd.read_file(f"zip://{zip_path}!conjunto_de_datos/{name}.shp")


# ---------------------------------------------------------------- simplificación
def simplify(gdf: gpd.GeoDataFrame, tol_m: float) -> gpd.GeoDataFrame:
    """Simplifica como cobertura (bordes compartidos idénticos → sin huecos entre vecinos)
    en metros (EPSG:6372) y reproyecta a lat/lon."""
    g = gdf.copy()
    g["geometry"] = shapely.coverage_simplify(g.geometry.values, tolerance=tol_m)
    return g.to_crs(4326)


def write_geojson(gdf, path: Path, precision: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(path, driver="GeoJSON", COORDINATE_PRECISION=precision)
    return path.stat().st_size / 1024


def main():
    zip_path = Path(sys.argv[1]) if len(sys.argv) > 1 else download_integrated()
    ent = read_layer(zip_path, "00ent")
    mun = read_layer(zip_path, "00mun")
    print(f"✔ {len(ent)} estados, {len(mun)} municipios (CRS EPSG:{ent.crs.to_epsg()})")

    # --- Web: estados
    e_web = simplify(ent[["CVE_ENT", "NOMGEO", "geometry"]], TOL_ESTADOS_M)
    kb = write_geojson(e_web, WEB_DIR / "estados.geojson", 3)
    print(f"✔ estados.geojson: {kb:.0f} KB")

    # --- Web: municipios por entidad
    for cve, g in mun.groupby("CVE_ENT"):
        write_geojson(simplify(g[["CVEGEO", "CVE_ENT", "CVE_MUN", "NOMGEO", "geometry"]], TOL_MUNICIPIOS_M),
                      WEB_DIR / "municipios" / f"{cve}.geojson", 4)
    print(f"✔ municipios/*.geojson: {mun.CVE_ENT.nunique()} archivos")

    # --- Catálogo para menús
    cat = {
        "fuente": SOURCE,
        "estados": [
            {"cve": r.CVE_ENT, "nombre": r.NOMGEO,
             "municipios": [{"cve": m.CVEGEO, "nombre": m.NOMGEO}
                            for m in mun[mun.CVE_ENT == r.CVE_ENT].sort_values("NOMGEO").itertuples()]}
            for r in ent.sort_values("NOMGEO").itertuples()
        ],
    }
    (WEB_DIR / "catalogo.json").write_text(json.dumps(cat, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print("✔ catalogo.json")

    # --- Cálculo zonal: geometrías más finas (50 m), una sola capa con nivel
    e_z = simplify(ent, TOL_ZONAL_M).assign(nivel="estado", CVEGEO=lambda d: d.CVE_ENT)
    m_z = simplify(mun, TOL_ZONAL_M).assign(nivel="municipio")
    zonal = pd.concat([e_z, m_z], ignore_index=True)[["CVEGEO", "CVE_ENT", "NOMGEO", "nivel", "geometry"]]
    config.ZONAL_GPKG.parent.mkdir(parents=True, exist_ok=True)
    zonal.to_file(config.ZONAL_GPKG, driver="GPKG")
    print(f"✔ {config.ZONAL_GPKG.name}: {len(zonal)} polígonos ({config.ZONAL_GPKG.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

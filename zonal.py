"""Estadística zonal ponderada por FRACCIÓN DE CELDA (la resolución de IMERG no cambia).

Cada celda de 0.1° pesa  w = (fracción de la celda dentro del polígono) × cos(lat)
(cos(lat) ∝ área real de la celda). El promedio de un día es  Σ w·P / Σ w  sobre las
celdas con dato ese día. Es el mismo criterio que usan exactextract / zonal stats.

    python zonal.py      # (re)calcula data/processed/zonal_series.nc sin reprocesar todo

Lo usan:
  - process.py  -> precalcula series de los 32 estados, 2,478 municipios y México ("00")
  - app/main.py -> áreas dibujadas por el usuario (polygon_weights al momento)
"""
import sys

import numpy as np
import shapely
import xarray as xr
from scipy import sparse

import config

KM_PER_DEG = 111.195  # km por grado de arco (R = 6371 km)


def grid_res(lat: np.ndarray) -> float:
    return float(np.round(np.median(np.diff(lat)), 6))


def polygon_weights(geom, lat: np.ndarray, lon: np.ndarray, res: float | None = None):
    """Índices (i_lat, j_lon) de las celdas que toca `geom` y la fracción cubierta (0–1].
    `lat`/`lon` son los centros de celda en orden ascendente."""
    res = res or grid_res(lat)
    h = res / 2
    w, s, e, n = geom.bounds
    i0, i1 = np.searchsorted(lat, s - h), np.searchsorted(lat, n + h)
    j0, j1 = np.searchsorted(lon, w - h), np.searchsorted(lon, e + h)
    if i0 >= i1 or j0 >= j1:
        return np.array([], int), np.array([], int), np.array([], float)
    jj, ii = np.meshgrid(np.arange(j0, j1), np.arange(i0, i1))
    ii, jj = ii.ravel(), jj.ravel()
    boxes = shapely.box(lon[jj] - h, lat[ii] - h, lon[jj] + h, lat[ii] + h)

    shapely.prepare(geom)
    hit = shapely.intersects(geom, boxes)
    ii, jj, boxes = ii[hit], jj[hit], boxes[hit]
    frac = np.ones(len(boxes))
    edge = ~shapely.contains_properly(geom, boxes)          # solo las del borde requieren intersección
    if edge.any():
        frac[edge] = shapely.area(shapely.intersection(boxes[edge], geom)) / (res * res)
    keep = frac > 1e-6
    return ii[keep], jj[keep], np.clip(frac[keep], 0, 1)


def weight_matrix(geoms, lat, lon) -> sparse.csr_matrix:
    """Matriz dispersa (n_regiones × n_celdas) con w = fracción × cos(lat)."""
    res, nlon = grid_res(lat), len(lon)
    coslat = np.cos(np.deg2rad(lat))
    rows, cols, vals = [], [], []
    for r, g in enumerate(geoms):
        ii, jj, frac = polygon_weights(g, lat, lon, res)
        rows.append(np.full(len(ii), r)); cols.append(ii * nlon + jj); vals.append(frac * coslat[ii])
    return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(len(geoms), len(lat) * nlon))


def compute_zonal_series() -> None:
    import geopandas as gpd

    files = sorted(config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"))
    if not files:
        sys.exit("No hay NetCDF procesados; corre process.py primero.")
    if not config.ZONAL_GPKG.exists():
        sys.exit(f"Falta {config.ZONAL_GPKG.name}; corre make_boundaries.py primero.")

    regions = gpd.read_file(config.ZONAL_GPKG)
    ds = xr.open_mfdataset(files, combine="by_coords")
    lat, lon = ds.lat.values, ds.lon.values
    res = grid_res(lat)

    print(f"Calculando pesos por fracción de celda para {len(regions)} polígonos…")
    W = weight_matrix(regions.geometry.values, lat, lon)
    # México completo = suma de los estados (no se traslapan)
    W = sparse.vstack([W[(regions.nivel == "estado").values].sum(axis=0), W]).tocsr()
    codes = np.concatenate([["00"], regions.CVEGEO.values])
    names = np.concatenate([["México"], regions.NOMGEO.values])
    levels = np.concatenate([["pais"], regions.nivel.values])
    ents = np.concatenate([["00"], regions.CVE_ENT.values])

    out, has_data = [], np.zeros(W.shape[1], bool)
    for year, dsy in ds.groupby("time.year"):
        P = dsy["precipitation"].values.reshape(dsy.sizes["time"], -1)   # (t, celdas)
        valid = ~np.isnan(P)
        has_data |= valid.any(0)
        num = (W @ np.nan_to_num(P).T)                                     # (regiones, t)
        den = (W @ valid.T.astype(np.float32))
        with np.errstate(invalid="ignore", divide="ignore"):
            out.append(np.where(den > 0, num / den, np.nan).astype(np.float32))
        print(f"  {year}: {dsy.sizes['time']} días")
    series = np.concatenate(out, axis=1)

    cell_km2 = (res * KM_PER_DEG) ** 2   # × cos(lat) ya incluido en W
    area_km2 = np.asarray(W @ has_data.astype(float)).ravel() * cell_km2

    zs = xr.Dataset(
        {"precipitation": (("region", "time"), series,
                           {"long_name": "Precipitación media diaria ponderada por fracción de celda", "units": "mm/day"})},
        coords={"region": codes, "time": ds.time.values, "nombre": ("region", names),
                "nivel": ("region", levels), "cve_ent": ("region", ents),
                "area_km2": ("region", np.round(area_km2, 1))},
        attrs={"method": "fracción de celda × cos(lat); renormalizado por celdas con dato",
               "boundaries": "INEGI, Marco Geoestadístico, Encuesta Intercensal 2025"},
    )
    enc = {"precipitation": {"zlib": True, "complevel": 4, "chunksizes": (min(256, len(codes)), series.shape[1])}}
    config.ZONAL_SERIES.unlink(missing_ok=True)
    zs.to_netcdf(config.ZONAL_SERIES, encoding=enc)
    n_empty = int(np.isnan(series).all(1).sum())
    print(f"✔ {config.ZONAL_SERIES.name}: {len(codes)} regiones × {series.shape[1]} días "
          f"({config.ZONAL_SERIES.stat().st_size / 1e6:.1f} MB); sin datos: {n_empty}")


if __name__ == "__main__":
    compute_zonal_series()

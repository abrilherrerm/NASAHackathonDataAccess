"""Paso 3: API FastAPI con series de tiempo interactivas.

    uvicorn app.main:app --reload
    -> http://127.0.0.1:8000        (mapa + gráfica interactiva)
    -> http://127.0.0.1:8000/docs   (Swagger)
"""
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import xarray as xr
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
import shapely
from shapely.geometry import shape

sys.path.append(str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402
from zonal import KM_PER_DEG, grid_res, polygon_weights  # noqa: E402

app = FastAPI(title="IMERG México – Precipitación diaria", version="1.0")
app.add_middleware(GZipMiddleware, minimum_size=1024)   # GeoJSON y JSON grandes viajan comprimidos
app.mount("/static", StaticFiles(directory=Path(__file__).parent / "static"), name="static")
Freq = Literal["D", "W", "MS", "YS"]
FREQ_LABEL = {"D": "diaria", "W": "semanal", "MS": "mensual", "YS": "anual"}


def _signature(paths) -> tuple:
    return tuple((str(p), p.stat().st_mtime) for p in paths if p.exists())


@lru_cache(maxsize=1)
def _open_daily(sig: tuple) -> xr.Dataset:
    # Abrir perezosamente respetando los chunks guardados en disco
    return xr.open_mfdataset([p for p, _ in sig], combine="by_coords", chunks={})


def dataset() -> xr.Dataset:
    """Se vuelve a abrir solo si process.py reescribió los archivos (no hace falta reiniciar)."""
    files = sorted(config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"))
    if not files:
        raise HTTPException(503, "No hay datos procesados. Corre download.py y process.py.")
    return _open_daily(_signature(files))


@lru_cache(maxsize=1)
def _open_zonal(sig: tuple) -> xr.Dataset:
    return xr.open_dataset(sig[0][0]).load()   # ~18 MB con 5 años: cabe en memoria


def zonal_dataset() -> xr.Dataset:
    if not config.ZONAL_SERIES.exists():
        raise HTTPException(503, "No hay series por entidad. Corre make_boundaries.py y zonal.py.")
    return _open_zonal(_signature([config.ZONAL_SERIES]))


@lru_cache(maxsize=1)
def _valid_cells(sig: tuple) -> np.ndarray:
    return dataset()["precipitation"].isel(time=slice(0, 31)).notnull().any("time").values


def _nearest_cell(lat: float, lon: float, max_km: float = 30) -> tuple[int, int]:
    """Índices de la celda CON DATO más cercana; 404 si no hay ninguna a menos de `max_km`."""
    pr = dataset()["precipitation"]
    valid = _valid_cells(_signature(sorted(config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"))))
    LA, LO = np.meshgrid(pr.lat.values, pr.lon.values, indexing="ij")
    d2 = np.where(valid, (LA - lat) ** 2 + ((LO - lon) * np.cos(np.deg2rad(lat))) ** 2, np.inf)
    i, j = np.unravel_index(np.argmin(d2), d2.shape)
    if np.sqrt(d2[i, j]) * KM_PER_DEG > max_km:
        raise HTTPException(404, "Punto fuera de México o sin datos.")
    return int(i), int(j)


def _slice(da: xr.DataArray, start: str | None, end: str | None) -> xr.DataArray:
    return da.sel(time=slice(start, end))


def late_start() -> pd.Timestamp | None:
    """Primer día con dato preliminar (IMERG Late), o None si todo es Final."""
    ds = dataset()
    if "is_late" not in ds:
        return None
    late = ds["is_late"].values == 1
    return pd.Timestamp(ds.time.values[late.argmax()]) if late.any() else None


def _late_from(s: pd.Series, freq: str) -> str | None:
    lf = late_start()
    if lf is None or s.empty:
        return None
    last = s.index.max()
    period_end = last if freq in ("D", "W") else last + pd.tseries.frequencies.to_offset(freq) - pd.Timedelta(days=1)
    return str(lf.date()) if period_end >= lf else None


def _series_payload(s: pd.Series, freq: str, **meta) -> dict:
    s = s.astype("float64").round(2)
    return {
        **meta,
        "late_from": _late_from(s, freq),
        "freq": freq,
        "units": "mm/día" if freq == "D" else "mm (acumulado)",
        "dates": s.index.strftime("%Y-%m-%d").tolist(),
        "values": [None if np.isnan(v) else float(v) for v in s.values],
        "stats": {
            "total_mm": round(float(np.nansum(s.values)), 1),
            "max": None if s.isna().all() else round(float(s.max()), 2),
            "max_date": None if s.isna().all() else s.idxmax().strftime("%Y-%m-%d"),
            "mean": None if s.isna().all() else round(float(s.mean()), 2),
        },
    }


def _resample(da: xr.DataArray, freq: str) -> pd.Series:
    s = da.to_series()
    # Diario: días sin datos quedan como null (huecos visibles, no líneas que los atraviesan)
    return s.asfreq("D") if freq == "D" else s.resample(freq).sum(min_count=1)


@app.get("/api/info")
def info():
    ds = dataset()
    t = ds.time.values
    return {
        "start": str(pd.Timestamp(t[0]).date()),
        "end": str(pd.Timestamp(t[-1]).date()),
        "n_days": int(t.size),
        "late_from": str(late_start().date()) if late_start() is not None else None,
        "bbox": config.MEXICO_BBOX,
        "resolution_deg": round(float(abs(ds.lat[1] - ds.lat[0])), 3),
    }


@app.get("/api/timeseries/point")
def point_series(
    lat: float = Query(..., ge=-90, le=90),
    lon: float = Query(..., ge=-180, le=180),
    start: str | None = None,
    end: str | None = None,
    freq: Freq = "D",
):
    """Serie en la celda IMERG (0.1°) con dato más cercana al punto (en la costa la exacta puede ser mar)."""
    i, j = _nearest_cell(lat, lon)
    da = _slice(dataset()["precipitation"].isel(lat=i, lon=j), start, end).load()
    if da.isnull().all():
        raise HTTPException(404, "Punto fuera de México o sin datos.")
    return _series_payload(
        _resample(da, freq), freq,
        lat=round(float(da.lat), 3), lon=round(float(da.lon), 3),
        label=f"Punto ({float(da.lat):.2f}, {float(da.lon):.2f})",
    )


@app.get("/api/timeseries/region")
def region_series(
    west: float = config.MEXICO_BBOX["west"], south: float = config.MEXICO_BBOX["south"],
    east: float = config.MEXICO_BBOX["east"], north: float = config.MEXICO_BBOX["north"],
    start: str | None = None, end: str | None = None, freq: Freq = "D",
):
    """Promedio espacial ponderado por cos(lat) dentro de un rectángulo (por defecto todo México)."""
    da = dataset()["precipitation"].sel(lat=slice(south, north), lon=slice(west, east))
    da = _slice(da, start, end)
    if da.size == 0:
        raise HTTPException(404, "Región vacía.")
    w = np.cos(np.deg2rad(da.lat))
    mean = da.weighted(w).mean(("lat", "lon")).load()
    return _series_payload(_resample(mean, freq), freq, label="Promedio regional",
                           bbox={"west": west, "south": south, "east": east, "north": north})


class PolygonRequest(BaseModel):
    geometry: dict[str, Any]          # GeoJSON: Polygon, MultiPolygon, Feature o FeatureCollection
    name: str | None = None
    start: str | None = None
    end: str | None = None
    freq: Freq = "D"


def _to_geometry(gj: dict):
    if gj.get("type") == "FeatureCollection":
        return shapely.union_all([shape(f["geometry"]) for f in gj["features"]])
    if gj.get("type") == "Feature":
        return shape(gj["geometry"])
    return shape(gj)


@app.post("/api/timeseries/polygon")
def polygon_series(req: PolygonRequest):
    """Precipitación promedio (ponderada por área) dentro de un polígono dibujado o de un estado.
    Además devuelve el área con datos (km²) y el volumen precipitado (hm³) en el periodo."""
    try:
        geom = _to_geometry(req.geometry)
    except Exception as e:
        raise HTTPException(422, f"GeoJSON inválido: {e}")
    if geom.is_empty or geom.area == 0:
        raise HTTPException(422, "El polígono está vacío.")
    if not geom.is_valid:
        geom = shapely.make_valid(geom)

    full = dataset()["precipitation"]
    ii, jj, frac = polygon_weights(geom, full.lat.values, full.lon.values)
    if len(ii) == 0:
        raise HTTPException(404, "El polígono está fuera del área de datos (México).")

    # Recorte al rectángulo de celdas tocadas + pesos = fracción × cos(lat)
    i0, i1, j0, j1 = ii.min(), ii.max() + 1, jj.min(), jj.max() + 1
    da = _slice(full.isel(lat=slice(i0, i1), lon=slice(j0, j1)), req.start, req.end)
    coslat = np.cos(np.deg2rad(da.lat.values))
    wgrid = np.zeros((i1 - i0, j1 - j0))
    wgrid[ii - i0, jj - j0] = frac * coslat[ii - i0]
    weights = xr.DataArray(wgrid, coords={"lat": da.lat, "lon": da.lon}, dims=("lat", "lon"))

    has_data = da.notnull().any("time").load()
    w_valid = float((weights * has_data).sum())
    if w_valid == 0:
        raise HTTPException(404, "El área no tiene datos (¿está sobre el mar o fuera de México?).")
    mean = da.weighted(weights).mean(("lat", "lon")).load()
    area_km2 = w_valid * (grid_res(full.lat.values) * KM_PER_DEG) ** 2
    n_cells = int(((weights > 0) & has_data).sum())

    payload = _series_payload(_resample(mean, req.freq), req.freq,
                              label=req.name or "Área seleccionada", n_cells=n_cells,
                              area_km2=round(area_km2, 1))
    total = payload["stats"]["total_mm"]
    payload["stats"]["volume_hm3"] = round(total * area_km2 / 1000, 1)   # mm·km² = 1,000 m³
    return payload


@app.get("/api/timeseries/entity/{cvegeo}")
def entity_series(cvegeo: str, start: str | None = None, end: str | None = None, freq: Freq = "D"):
    """Serie precalculada de un estado (CVE_ENT, p. ej. '20'), municipio (CVEGEO, p. ej. '20067')
    o de México completo ('00'). Promedio ponderado por fracción de celda."""
    zs = zonal_dataset()
    if cvegeo not in zs.region.values:
        raise HTTPException(404, f"Clave {cvegeo!r} no encontrada.")
    r = zs.sel(region=cvegeo)
    da = _slice(r["precipitation"], start, end)
    area = float(r.area_km2)
    nivel = str(r.nivel.values)
    nombre = str(r.nombre.values)
    if nivel == "municipio":
        ent = str(zs.sel(region=str(r.cve_ent.values)).nombre.values)
        label = f"{nombre}, {ent}"
    else:
        label = nombre
    payload = _series_payload(_resample(da, freq), freq, label=label, cvegeo=cvegeo, nivel=nivel,
                              area_km2=round(area, 1))
    payload["stats"]["volume_hm3"] = round(payload["stats"]["total_mm"] * area / 1000, 1)
    return payload


@app.get("/api/map")
def precip_map(start: str, end: str | None = None, step: int = Query(2, ge=1, le=10)):
    """Acumulado en [start, end] para pintar el mapa (submuestreado cada `step` celdas)."""
    da = _slice(dataset()["precipitation"], start, end or start)
    if da.sizes["time"] == 0:
        raise HTTPException(404, f"Sin datos procesados entre {start} y {end or start}. "
                                 f"Si ya se descargaron, corre process.py.")
    acc = da.sum("time", min_count=1).isel(lat=slice(None, None, step), lon=slice(None, None, step)).load()
    z = np.where(np.isnan(acc.values), None, np.round(acc.values, 1))
    return {"lat": acc.lat.values.round(3).tolist(), "lon": acc.lon.values.round(3).tolist(),
            "z": z.tolist(), "start": start, "end": end or start}


# ---------------------------------------------------------------------------
# Medidas para análisis autoguiado (measures.py)
# ---------------------------------------------------------------------------
import json  # noqa: E402

import measures  # noqa: E402
from analysis import _json_default  # noqa: E402


def _daily_for_place(place: dict) -> tuple[pd.Series, str]:
    """Serie diaria (mm/día) para un lugar, con el calendario completo (días faltantes = NaN)."""
    if "cve" in place:
        zs = zonal_dataset()
        if place["cve"] not in zs.region.values:
            raise HTTPException(404, f"Clave {place['cve']!r} no encontrada.")
        r = zs.sel(region=place["cve"])
        s, label = r["precipitation"].to_series(), str(r.nombre.values)
    elif "point" in place:
        i, j = _nearest_cell(*place["point"])
        da = dataset()["precipitation"].isel(lat=i, lon=j).load()
        s, label = da.to_series(), f"Celda ({float(da.lat):.2f}, {float(da.lon):.2f})"
    elif "geometry" in place:
        geom = _to_geometry(place["geometry"])
        full = dataset()["precipitation"]
        ii, jj, frac = polygon_weights(geom, full.lat.values, full.lon.values)
        if len(ii) == 0:
            raise HTTPException(404, "El polígono está fuera del área de datos (México).")
        i0, i1, j0, j1 = ii.min(), ii.max() + 1, jj.min(), jj.max() + 1
        da = full.isel(lat=slice(i0, i1), lon=slice(j0, j1))
        wgrid = np.zeros((i1 - i0, j1 - j0))
        wgrid[ii - i0, jj - j0] = frac * np.cos(np.deg2rad(da.lat.values))[ii - i0]
        w = xr.DataArray(wgrid, coords={"lat": da.lat, "lon": da.lon}, dims=("lat", "lon"))
        s, label = da.weighted(w).mean(("lat", "lon")).load().to_series(), place.get("name") or "Área"
    else:
        raise HTTPException(422, "Lugar inválido: usa cve, point o geometry.")
    s = s.astype("float64")
    return s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="D")), label


measures.set_daily_source(_daily_for_place)


class MeasureRequest(BaseModel):
    measure: str
    place: dict[str, Any]
    months: list[int] = list(range(1, 13))
    period: list[int] | None = None          # [año_inicio, año_fin]
    baseline: list[int] | None = None
    threshold: float = 1.0                   # wet_days
    dias: int = 7                            # dry_spell
    scale: int = 3                           # drought
    a: list[int] | None = None               # compare_periods
    b: list[int] | None = None
    place2: dict[str, Any] | None = None     # compare_places


@app.get("/api/measures")
def list_measures():
    """Catálogo de medidas con su disponibilidad según los datos descargados."""
    out = []
    for m in measures.CATALOG.values():
        ok, why = measures.available(m)
        out.append({**measures.asdict(m), "disponible": ok, "motivo": why})
    return out


@app.post("/api/measure")
def run_measure(req: MeasureRequest):
    m = measures.CATALOG.get(req.measure)
    if not m:
        raise HTTPException(404, f"Medida {req.measure!r} no existe.")
    ok, why = measures.available(m)
    if not ok:
        raise HTTPException(503, why)
    if not req.months or any(not 1 <= x <= 12 for x in req.months):
        raise HTTPException(422, "Meses inválidos (1–12).")
    params = req.model_dump(exclude={"measure", "place"}, exclude_none=True)
    for k in ("period", "baseline", "a", "b"):
        if k in params:
            params[k] = tuple(params[k])
    try:
        res = measures.compute(req.measure, req.place, **params)
    except HTTPException:
        raise
    except (ValueError, KeyError) as e:
        raise HTTPException(422, str(e))
    return _clean(json.loads(json.dumps(res, default=_json_default)))


def _clean(o):
    """NaN/inf -> null (JSON estricto) y claves numéricas (años) como texto."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(Path(__file__).parent / "static" / "index.html")

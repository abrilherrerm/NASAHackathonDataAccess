"""Paso 3: API FastAPI con series de tiempo interactivas.

    uvicorn app.main:app --reload
    -> http://127.0.0.1:8000        (mapa + gráfica interactiva)
    -> http://127.0.0.1:8000/docs   (Swagger)
"""
import sys
from functools import lru_cache
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import xarray as xr
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

sys.path.append(str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402

app = FastAPI(title="IMERG México – Precipitación diaria", version="1.0")
app.mount("/static", StaticFiles(directory=Path(__file__).parent / "static"), name="static")
Freq = Literal["D", "W", "MS", "YS"]
FREQ_LABEL = {"D": "diaria", "W": "semanal", "MS": "mensual", "YS": "anual"}


@lru_cache(maxsize=1)
def dataset() -> xr.Dataset:
    files = sorted(config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"))
    if not files:
        raise HTTPException(503, "No hay datos procesados. Corre download.py y process.py.")
    # Abrir perezosamente respetando los chunks guardados en disco
    return xr.open_mfdataset(files, combine="by_coords", chunks={})


def _slice(da: xr.DataArray, start: str | None, end: str | None) -> xr.DataArray:
    return da.sel(time=slice(start, end))


def _series_payload(s: pd.Series, freq: str, **meta) -> dict:
    s = s.astype("float64").round(2)
    return {
        **meta,
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
    return s if freq == "D" else s.resample(freq).sum(min_count=1)


@app.get("/api/info")
def info():
    ds = dataset()
    t = ds.time.values
    return {
        "start": str(pd.Timestamp(t[0]).date()),
        "end": str(pd.Timestamp(t[-1]).date()),
        "n_days": int(t.size),
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
    """Serie en la celda IMERG (0.1°) más cercana al punto."""
    da = dataset()["precipitation"].sel(lat=lat, lon=lon, method="nearest")
    da = _slice(da, start, end).load()
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


@app.get("/api/map")
def precip_map(start: str, end: str | None = None, step: int = Query(2, ge=1, le=10)):
    """Acumulado en [start, end] para pintar el mapa (submuestreado cada `step` celdas)."""
    da = _slice(dataset()["precipitation"], start, end or start)
    if da.sizes["time"] == 0:
        raise HTTPException(404, "Sin datos en esas fechas.")
    acc = da.sum("time", min_count=1).isel(lat=slice(None, None, step), lon=slice(None, None, step)).load()
    z = np.where(np.isnan(acc.values), None, np.round(acc.values, 1))
    return {"lat": acc.lat.values.round(3).tolist(), "lon": acc.lon.values.round(3).tolist(),
            "z": z.tolist(), "start": start, "end": end or start}


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(Path(__file__).parent / "static" / "index.html")

"""IMERG Final V07 MENSUAL (GPM_3IMERGM) 1998 → hoy, recortado a México: descarga + proceso.

Complementa la serie diaria (2020→): da ~27 años de totales mensuales para climatologías de
referencia, SPI, anomalías y tendencias. Las métricas diarias (días con lluvia, rachas secas,
Rx1day) siguen necesitando el producto diario.

    python download_monthly.py               # descarga lo que falte (1 job de Harmony por año) y procesa
    python download_monthly.py --process     # sólo reprocesa lo ya descargado

Salidas:
    data/processed/imerg_mexico_monthly.nc   precipitación (mm/mes), error aleatorio (mm/mes),
                                             índice de calidad (pluviómetros equivalentes) y peso de pluviómetros (%)
    data/processed/zonal_monthly.nc          las mismas series por México / estado / municipio (INEGI)

OJO: en el archivo original `precipitation` y `randomError` vienen en mm/h (tasa media del mes);
aquí se convierten a mm/mes multiplicando por 24 × días del mes.
"""
import argparse
import re
import sys

import numpy as np
import pandas as pd
import xarray as xr
from harmony import Collection, Request

import config
from download import last_available, login, mexico_geojson
from process import mexico_mask

VARS = {"precipitation": "precipitation", "randomError": "random_error",
        "precipitationQualityIndex": "quality_index", "gaugeRelativeWeighting": "gauge_weight"}
DATE_RE = re.compile(r"3IMERG\.(\d{8})-")


def download(start: str, end: str) -> None:
    harmony = login()
    shape = mexico_geojson()
    collection = Collection(id=config.MONTHLY_CONCEPT_ID)
    for year in range(pd.Timestamp(start).year, pd.Timestamp(end).year + 1):
        t0 = max(pd.Timestamp(f"{year}-01-01"), pd.Timestamp(start))
        t1 = min(pd.Timestamp(f"{year}-12-31 23:59:59"), pd.Timestamp(end) + pd.Timedelta(days=1, seconds=-1))
        out_dir = config.RAW_MONTHLY_DIR / str(year)
        expected = len(pd.period_range(t0, t1, freq="M"))
        present = {DATE_RE.search(f.name).group(1) for f in out_dir.glob("*.nc*") if DATE_RE.search(f.name)}
        if len(present) >= expected:
            print(f"↷ {year} completo, se omite")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        for partial in out_dir.glob("*.nc*"):
            partial.unlink()
        request = Request(collection=collection, shape=str(shape),
                          temporal={"start": t0.to_pydatetime(), "stop": t1.to_pydatetime()},
                          variables=[f"Grid/{v}" for v in VARS])
        if not request.is_valid():
            sys.exit(f"Petición inválida: {request.error_messages()}")
        job_id = harmony.submit(request)
        print(f"→ {year}: job Harmony {job_id}")
        harmony.wait_for_processing(job_id, show_progress=True)
        files = [f.result() for f in harmony.download_all(job_id, directory=str(out_dir), overwrite=True)]
        print(f"✔ {year}: {len(files)} archivos")


def process() -> None:
    files = sorted(config.RAW_MONTHLY_DIR.rglob("*.nc*"))
    if not files:
        sys.exit(f"No hay archivos en {config.RAW_MONTHLY_DIR}.")
    by_month = {}
    for f in files:   # si un mes se bajó dos veces, gana el más reciente
        m = DATE_RE.search(f.name)
        if m and (m.group(1) not in by_month or f.stat().st_mtime > by_month[m.group(1)].stat().st_mtime):
            by_month[m.group(1)] = f
    parts = []
    for date, f in sorted(by_month.items()):
        # Las variables viven en el grupo HDF5 "Grid"; el tiempo viene en cftime -> se toma del nombre
        with xr.open_dataset(f, engine="h5netcdf", group="Grid", decode_times=False) as g:
            keep = [v for v in VARS if v in g]
            d = g[keep].rename({k: VARS[k] for k in keep}).load()
        d = d.assign_coords(time=[pd.Timestamp(date)])
        parts.append(d.transpose("time", "lat", "lon"))
    ds = xr.concat(parts, dim="time").sortby("lat")
    hours = xr.DataArray(ds.time.dt.days_in_month.values * 24.0, coords={"time": ds.time}, dims="time")
    for v in ("precipitation", "random_error"):
        if v in ds:
            ds[v] = (ds[v] * hours).astype("float32")
            ds[v].attrs = {"units": "mm/month"}
    ds["precipitation"].attrs["long_name"] = "Precipitación total mensual (IMERG Final V07, satélite + pluviómetros)"
    if "random_error" in ds:
        ds["random_error"].attrs["long_name"] = "Error aleatorio estimado de la precipitación mensual"
    if "quality_index" in ds:
        ds["quality_index"].attrs = {"long_name": "Índice de calidad", "units": "pluviómetros equivalentes por caja de 2.5°"}
    if "gauge_weight" in ds:
        ds["gauge_weight"].attrs = {"long_name": "Peso relativo de pluviómetros", "units": "%"}
    ds = ds.map(lambda da: mexico_mask(da) if {"lat", "lon"} <= set(da.dims) else da)
    ds.attrs = {"title": "Precipitación mensual IMERG V07 recortada a México",
                "source": f"GPM_3IMERGM_07 ({config.MONTHLY_CONCEPT_ID}) vía NASA Harmony"}
    enc = {v: {"zlib": True, "complevel": 4, "dtype": "float32", "_FillValue": np.float32(-9999.9)} for v in ds.data_vars}
    config.MONTHLY_FILE.unlink(missing_ok=True)
    ds.to_netcdf(config.MONTHLY_FILE, encoding=enc)
    t = ds.time.to_index()
    print(f"✔ {config.MONTHLY_FILE.name}: {t[0]:%Y-%m} → {t[-1]:%Y-%m} ({len(t)} meses), variables {list(ds.data_vars)}")
    zonal_monthly(ds)


def _regions():
    """Polígonos INEGI: el GPKG de cálculo (50 m) si existe; si no, los GeoJSON del frontend (250 m),
    suficientes para celdas de ~11 km."""
    import geopandas as gpd
    if config.ZONAL_GPKG.exists():
        return gpd.read_file(config.ZONAL_GPKG)
    bdir = config.BASE_DIR / "app" / "static" / "boundaries"
    est = gpd.read_file(bdir / "estados.geojson").assign(nivel="estado")
    est["CVEGEO"] = est["CVE_ENT"]
    mun = pd.concat([gpd.read_file(f) for f in sorted((bdir / "municipios").glob("*.geojson"))], ignore_index=True)
    mun = mun.assign(nivel="municipio")
    mun["CVE_ENT"] = mun["CVEGEO"].str[:2]
    return pd.concat([est, mun], ignore_index=True)[["CVEGEO", "NOMGEO", "CVE_ENT", "nivel", "geometry"]]


def zonal_monthly(ds: xr.Dataset) -> None:
    from scipy import sparse
    from zonal import KM_PER_DEG, grid_res, weight_matrix
    regions = _regions()
    lat, lon = ds.lat.values, ds.lon.values
    print(f"Pesos por fracción de celda para {len(regions)} polígonos…")
    W = weight_matrix(regions.geometry.values, lat, lon)
    W = sparse.vstack([W[(regions.nivel == "estado").values].sum(axis=0), W]).tocsr()
    data = {}
    for v in ("precipitation", "random_error"):
        if v not in ds:
            continue
        P = ds[v].values.reshape(ds.sizes["time"], -1)
        valid = ~np.isnan(P)
        num, den = W @ np.nan_to_num(P).T, W @ valid.T.astype(np.float32)
        with np.errstate(invalid="ignore", divide="ignore"):
            data[v] = (("region", "time"), np.where(den > 0, num / den, np.nan).astype(np.float32), {"units": "mm/month"})
    has = ~np.isnan(ds["precipitation"].values).all(0).ravel()
    area = np.asarray(W @ has.astype(float)).ravel() * (grid_res(lat) * KM_PER_DEG) ** 2
    zs = xr.Dataset(data, coords={
        "region": np.concatenate([["00"], regions.CVEGEO.values]), "time": ds.time.values,
        "nombre": ("region", np.concatenate([["México"], regions.NOMGEO.values])),
        "nivel": ("region", np.concatenate([["pais"], regions.nivel.values])),
        "cve_ent": ("region", np.concatenate([["00"], regions.CVE_ENT.values])),
        "area_km2": ("region", np.round(area, 1))},
        attrs={"method": "fracción de celda × cos(lat)",
               "random_error": "promedio ponderado del error por celda = cota superior (supone errores "
                               "totalmente correlacionados entre celdas)"})
    config.ZONAL_MONTHLY.unlink(missing_ok=True)
    zs.to_netcdf(config.ZONAL_MONTHLY)
    print(f"✔ {config.ZONAL_MONTHLY.name}: {zs.sizes['region']} regiones × {zs.sizes['time']} meses")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--start", default=config.MONTHLY_START)
    p.add_argument("--end", default=config.END_DATE)
    p.add_argument("--process", action="store_true", help="sólo procesar lo descargado")
    args = p.parse_args()
    if not args.process:
        # END_DATE = None en config -> lo último publicado en CMR para el producto mensual
        end = args.end or str(last_available(config.MONTHLY_CONCEPT_ID).date())
        download(args.start, end)
    process()


if __name__ == "__main__":
    main()

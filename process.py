"""Paso 2: une los recortes diarios de Harmony y guarda NetCDF anuales con chunks.

Harmony ya recorta con el polígono de México (fuera queda como _FillValue); aquí además
se aplica la misma máscara localmente como respaldo.

Uso:
    python process.py
    python process.py --no-mask   # no reaplicar máscara local
"""
import argparse
import re
import sys

import numpy as np
import xarray as xr

import config


def open_raw() -> xr.DataArray:
    files = sorted(config.RAW_DIR.rglob("*.nc*"))
    if not files:
        sys.exit(f"No hay archivos en {config.RAW_DIR}. Corre download.py primero.")
    files = dedupe_by_date(files)
    ok = []
    for f in files:   # descarta archivos truncados (p. ej. si se cortó la descarga)
        try:
            with xr.open_dataset(f, engine="h5netcdf") as tmp:
                tmp[config.VARIABLE].shape
            ok.append(f)
        except Exception as e:
            print(f"⚠ se omite {f.name}: {e}")
    if not ok:
        sys.exit("✘ Ningún archivo se pudo leer (ver mensajes arriba).")
    print(f"Abriendo {len(ok)} de {len(files)} archivos…")
    ds = xr.open_mfdataset(ok, combine="by_coords", chunks={"time": 31}, engine="h5netcdf")
    pr = ds[config.VARIABLE].transpose("time", "lat", "lon")   # IMERG viene (time, lon, lat)
    pr["time"] = pr["time"].dt.floor("D")
    pr = pr.sortby("lat").sortby("time")
    pr.name = "precipitation"
    pr.attrs = {"long_name": "Precipitación total diaria (IMERG Final V07)", "units": "mm/day"}
    return pr


def dedupe_by_date(files):
    """Si un día se descargó dos veces (descarga reanudada), usa solo el archivo más reciente."""
    by_date = {}
    for f in files:
        m = re.search(r"3IMERG\.(\d{8})-", f.name)
        key = m.group(1) if m else f.name
        if key not in by_date or f.stat().st_mtime > by_date[key].stat().st_mtime:
            by_date[key] = f
    if len(by_date) < len(files):
        print(f"ℹ {len(files) - len(by_date)} archivos duplicados ignorados (mismo día)")
    return sorted(by_date.values())


def mexico_mask(pr: xr.DataArray) -> xr.DataArray:
    import regionmask
    countries = regionmask.defined_regions.natural_earth_v5_1_2.countries_50
    mask = countries.mask(pr.lon, pr.lat)
    return pr.where(mask == countries.map_keys("Mexico"))


def save_yearly(pr: xr.DataArray) -> None:
    config.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    for old in config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"):   # quita sintéticos/viejos
        old.unlink()
    ds = pr.to_dataset()
    ds.attrs.update({
        "title": "Precipitación diaria IMERG V07 recortada a México",
        "source": f"{config.SHORT_NAME}_{config.VERSION} ({config.CONCEPT_ID}) vía NASA Harmony",
    })
    for year, dsy in ds.groupby("time.year"):
        chunks = {k: min(v, dsy.sizes[k]) for k, v in config.NC_CHUNKS.items()}
        enc = {"precipitation": {
            "zlib": True, "complevel": 4, "dtype": "float32", "_FillValue": np.float32(-9999.9),
            "chunksizes": (chunks["time"], chunks["lat"], chunks["lon"]),
        }}
        out = config.PROCESSED_DIR / f"imerg_mexico_daily_{year}.nc"
        dsy.chunk(chunks).to_netcdf(out, encoding=enc, engine="netcdf4")
        print(f"✔ {out.name}: {dict(dsy.sizes)}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--no-mask", action="store_true")
    args = p.parse_args()
    pr = open_raw()
    save_yearly(pr if args.no_mask else mexico_mask(pr))

    # Series precalculadas por estado/municipio (fracción de celda)
    if config.ZONAL_GPKG.exists():
        from zonal import compute_zonal_series
        compute_zonal_series()
    else:
        print("ℹ Sin límites INEGI: corre make_boundaries.py y luego zonal.py para las series por entidad.")


if __name__ == "__main__":
    main()

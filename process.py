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
import pandas as pd
import xarray as xr

import config


LATE_DIR = config.RUNS["late"]["raw_dir"]


def _is_late(path) -> bool:
    return LATE_DIR in path.parents


def open_raw() -> xr.Dataset:
    """Une Final (data/raw) y Late (data/raw_late). Si un día existe en ambos, gana Final.
    Devuelve `precipitation` (time, lat, lon) y `is_late` (time): 1 = dato preliminar (Late)."""
    files = sorted(config.RAW_DIR.rglob("*.nc*")) + sorted(LATE_DIR.rglob("*.nc*"))
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
    n_late = sum(_is_late(f) for f in ok)
    print(f"Abriendo {len(ok)} de {len(files)} archivos ({len(ok) - n_late} Final, {n_late} Late)…")

    late_set = {str(f) for f in ok if _is_late(f)}

    def tag(ds):
        src = ds.encoding.get("source", "")
        flag = np.int8(1 if src in late_set else 0)
        return ds[[config.VARIABLE]].assign(is_late=("time", np.full(ds.sizes["time"], flag)))

    ds = xr.open_mfdataset(ok, combine="by_coords", chunks={"time": 31}, engine="h5netcdf", preprocess=tag)
    pr = ds[config.VARIABLE].transpose("time", "lat", "lon")   # IMERG viene (time, lon, lat)
    out = xr.Dataset({"precipitation": pr, "is_late": ds["is_late"]})
    out["time"] = out["time"].dt.floor("D")
    out = out.sortby("lat").sortby("time")
    out["precipitation"].attrs = {"long_name": "Precipitación total diaria (IMERG V07)", "units": "mm/day"}
    out["is_late"].attrs = {"long_name": "Fuente del día", "flag_values": "0, 1",
                            "flag_meanings": "final late", "comment": "1 = IMERG Late (preliminar)"}
    return out


def dedupe_by_date(files):
    """Un archivo por día: prefiere Final sobre Late y, a igualdad, el más reciente
    (descargas reanudadas pueden dejar el mismo día dos veces)."""
    by_date = {}
    rank = lambda f: (0 if _is_late(f) else 1, f.stat().st_mtime)
    for f in files:
        m = re.search(r"3IMERG\.(\d{8})-", f.name)
        key = m.group(1) if m else f.name
        if key not in by_date or rank(f) > rank(by_date[key]):
            by_date[key] = f
    if len(by_date) < len(files):
        print(f"ℹ {len(files) - len(by_date)} archivos ignorados (día repetido o ya cubierto por Final)")
    return sorted(by_date.values())


def mexico_mask(pr: xr.DataArray) -> xr.DataArray:
    import regionmask
    countries = regionmask.defined_regions.natural_earth_v5_1_2.countries_50
    mask = countries.mask(pr.lon, pr.lat)
    return pr.where(mask == countries.map_keys("Mexico"))


def save_yearly(ds: xr.Dataset) -> None:
    config.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    for old in config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"):   # quita sintéticos/viejos
        old.unlink()
    ds.attrs.update({
        "title": "Precipitación diaria IMERG V07 recortada a México",
        "source": "GPM_3IMERGDF_07 (Final) y GPM_3IMERGDL_07 (Late, solo donde no hay Final) vía NASA Harmony",
    })
    for year, dsy in ds.groupby("time.year"):
        chunks = {k: min(v, dsy.sizes[k]) for k, v in config.NC_CHUNKS.items()}
        enc = {"precipitation": {
            "zlib": True, "complevel": 4, "dtype": "float32", "_FillValue": np.float32(-9999.9),
            "chunksizes": (chunks["time"], chunks["lat"], chunks["lon"]),
        }, "is_late": {"dtype": "int8"}}
        out = config.PROCESSED_DIR / f"imerg_mexico_daily_{year}.nc"
        dsy.chunk(chunks).to_netcdf(out, encoding=enc, engine="netcdf4")
        print(f"✔ {out.name}: {dict(dsy.sizes)}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--no-mask", action="store_true")
    args = p.parse_args()
    ds = open_raw()
    if not args.no_mask:
        ds["precipitation"] = mexico_mask(ds["precipitation"])
    save_yearly(ds)

    # Series precalculadas por estado/municipio (fracción de celda)
    if config.ZONAL_GPKG.exists():
        from zonal import compute_zonal_series
        compute_zonal_series()
    else:
        print("ℹ Sin límites INEGI: corre make_boundaries.py y luego zonal.py para las series por entidad.")


if __name__ == "__main__":
    main()

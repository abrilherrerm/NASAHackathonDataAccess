"""Paso 1: login con earthaccess + recorte de México (polígono) con Harmony, por chunks mensuales.

Credenciales (NASA Earthdata Login), en este orden:
  1. Archivo .env con EARTHDATA_USERNAME / EARTHDATA_PASSWORD   (ver .env.example)
  2. ~/.netrc (se crea solo la primera vez que usas la opción 3)
  3. Si no hay ninguno, te las pide en la terminal y las guarda en ~/.netrc

Uso:
    python download.py                                      # 5 años completos
    python download.py --start 2025-01-01 --end 2025-01-31  # prueba de 1 mes
"""
import argparse
import os
import sys
from pathlib import Path

import earthaccess
import pandas as pd
from dotenv import load_dotenv
from harmony import Client, Collection, Request

import config


def login() -> Client:
    load_dotenv(config.BASE_DIR / ".env")
    user, pwd = os.getenv("EARTHDATA_USERNAME"), os.getenv("EARTHDATA_PASSWORD")

    if user and pwd:
        auth = earthaccess.login(strategy="environment")
        harmony = Client(auth=(user, pwd))
    else:
        has_netrc = (Path.home() / ".netrc").exists()
        auth = earthaccess.login(strategy="netrc" if has_netrc else "interactive", persist=True)
        harmony = Client()  # lee ~/.netrc
    if not auth or not auth.authenticated:
        sys.exit("✘ No se pudo autenticar con Earthdata Login (revisa usuario/contraseña).")
    print("✔ Autenticado en Earthdata Login")
    return harmony


def mexico_geojson() -> Path:
    """Polígono de México (Natural Earth 1:50m, simplificado) para el recorte en Harmony."""
    path = config.MEXICO_GEOJSON
    if not path.exists():
        import geopandas as gpd
        import regionmask
        countries = regionmask.defined_regions.natural_earth_v5_1_2.countries_50
        poly = countries[countries.map_keys("Mexico")].polygon.simplify(0.02)
        path.parent.mkdir(parents=True, exist_ok=True)
        gpd.GeoDataFrame({"name": ["Mexico"]}, geometry=[poly], crs="EPSG:4326").to_file(path, driver="GeoJSON")
        print(f"✔ Polígono de México guardado en {path}")
    return path


def month_chunks(start: str, end: str):
    """Divide el periodo en ventanas mensuales (chunks) para Harmony."""
    for m in pd.period_range(start, end, freq="M"):
        t0 = max(m.start_time, pd.Timestamp(start))
        t1 = min(m.end_time.floor("s"), pd.Timestamp(end) + pd.Timedelta(days=1, seconds=-1))
        yield t0.to_pydatetime(), t1.to_pydatetime()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--start", default=config.START_DATE)
    p.add_argument("--end", default=config.END_DATE)
    args = p.parse_args()

    harmony = login()
    shape = mexico_geojson()
    collection = Collection(id=config.CONCEPT_ID)

    for t0, t1 in month_chunks(args.start, args.end):
        out_dir = config.RAW_DIR / t0.strftime("%Y-%m")
        expected_days = (t1.date() - t0.date()).days + 1
        if out_dir.exists() and len(list(out_dir.glob("*.nc*"))) >= expected_days:
            print(f"↷ {t0:%Y-%m} completo, se omite")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)

        request = Request(
            collection=collection,
            shape=str(shape),               # recorte con el polígono de México
            temporal={"start": t0, "stop": t1},
            variables=[config.VARIABLE],
        )
        if not request.is_valid():
            sys.exit(f"Petición inválida: {request.error_messages()}")

        job_id = harmony.submit(request)
        print(f"→ {t0:%Y-%m}: job Harmony {job_id}")
        harmony.wait_for_processing(job_id, show_progress=True)
        files = [f.result() for f in harmony.download_all(job_id, directory=str(out_dir), overwrite=True)]
        print(f"✔ {t0:%Y-%m}: {len(files)} archivos")


if __name__ == "__main__":
    main()

"""Paso 1: login con earthaccess + recorte de México (polígono) con Harmony, por chunks mensuales.

Credenciales (NASA Earthdata Login), en este orden:
  1. Archivo .env con EARTHDATA_USERNAME / EARTHDATA_PASSWORD   (ver .env.example)
  2. ~/.netrc (se crea solo la primera vez que usas la opción 3)
  3. Si no hay ninguno, te las pide en la terminal y las guarda en ~/.netrc

Uso:
    python download.py                                      # 2014-06 → hoy (Final + Late)
    python download.py --start 2025-01-01 --end 2025-01-31  # prueba de 1 mes
    python download.py --run final                          # solo Final (sin datos preliminares)
"""
import argparse
import os
import re
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


def last_available(concept_id: str) -> pd.Timestamp:
    """Último día disponible: fin declarado de la colección o, si sigue activa, último granulo."""
    col = earthaccess.search_datasets(concept_id=concept_id)[0]["umm"]
    rng = col["TemporalExtents"][0]["RangeDateTimes"][0]
    if rng.get("EndingDateTime") and not col["TemporalExtents"][0].get("EndsAtPresentFlag"):
        return pd.Timestamp(rng["EndingDateTime"][:10])
    today = pd.Timestamp.today().normalize()
    grans = earthaccess.search_data(concept_id=concept_id,
                                    temporal=((today - pd.Timedelta(days=60)).strftime("%Y-%m-%d"),
                                              today.strftime("%Y-%m-%d")))
    if not grans:
        sys.exit(f"Sin granulos recientes para {concept_id}")
    return max(pd.Timestamp(g["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"][:10]) for g in grans)


def plan(start: str, end: str | None, run: str):
    """Tramos (run, inicio, fin): Final hasta su último día; Late después (si run='auto')."""
    final_end = last_available(config.RUNS["final"]["concept_id"])
    late_end = last_available(config.RUNS["late"]["concept_id"])
    end_ts = pd.Timestamp(end) if end else (late_end if run != "final" else final_end)
    start_ts = pd.Timestamp(start)
    print(f"ℹ Final disponible hasta {final_end.date()}, Late hasta {late_end.date()}")
    segs = []
    if run in ("auto", "final") and start_ts <= min(end_ts, final_end):
        segs.append(("final", start_ts, min(end_ts, final_end)))
    if run == "late":
        segs.append(("late", start_ts, min(end_ts, late_end)))
    elif run == "auto" and end_ts > final_end:
        segs.append(("late", max(start_ts, final_end + pd.Timedelta(days=1)), min(end_ts, late_end)))
    return segs


def download_segment(harmony, shape, run: str, start, end) -> None:
    cfg = config.RUNS[run]
    collection = Collection(id=cfg["concept_id"])
    for t0, t1 in month_chunks(start, end):
        out_dir = cfg["raw_dir"] / t0.strftime("%Y-%m")
        expected_days = (t1.date() - t0.date()).days + 1
        days_present = {re.search(r"3IMERG\.(\d{8})-", f.name).group(1)
                        for f in out_dir.glob("*.nc*") if re.search(r"3IMERG\.(\d{8})-", f.name)}
        if len(days_present) >= expected_days:
            print(f"↷ [{run}] {t0:%Y-%m} completo, se omite")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        for partial in out_dir.glob("*.nc*"):   # restos de un intento anterior incompleto
            partial.unlink()

        request = Request(
            collection=collection,
            shape=str(shape),               # recorte con el polígono de México
            temporal={"start": t0, "stop": t1},
            variables=[config.VARIABLE],
        )
        if not request.is_valid():
            sys.exit(f"Petición inválida: {request.error_messages()}")

        job_id = harmony.submit(request)
        print(f"→ [{run}] {t0:%Y-%m}: job Harmony {job_id}")
        harmony.wait_for_processing(job_id, show_progress=True)
        files = [f.result() for f in harmony.download_all(job_id, directory=str(out_dir), overwrite=True)]
        print(f"✔ [{run}] {t0:%Y-%m}: {len(files)} archivos")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--start", default=config.START_DATE)
    p.add_argument("--end", default=config.END_DATE, help="por defecto: lo más reciente disponible")
    p.add_argument("--run", choices=["auto", "final", "late"], default="auto",
                   help="auto = Final donde exista y Late después (por defecto)")
    args = p.parse_args()

    harmony = login()
    shape = mexico_geojson()
    for run, t0, t1 in plan(args.start, args.end, args.run):
        print(f"=== {run.upper()}: {t0.date()} → {t1.date()} ===")
        download_segment(harmony, shape, run, t0, t1)


if __name__ == "__main__":
    main()

"""¿Se puede integrar IMERG Late (preliminar) con Final en la estadística?

Late no tiene la calibración mensual con pluviómetros (GPCC) que sí tiene Final. Para saber si sus
datos (2025-10 → hoy) pueden sumarse a la serie Final, se compara Late contra Final en un periodo
donde existen ambos (por defecto oct-2023 → sep-2025, dos años completos):

  1. Sesgo de totales (mensual y por temporada) por mes calendario y por región climática.
  2. Detección de días con lluvia (≥ 1 mm): frecuencia, POD, FAR, CSI.
  3. Equivalencia (TOST) con margen ±10 %: "compatible" exige DEMOSTRAR que la diferencia es
     pequeña, no sólo "no encontrar diferencia" (eso sería confundir falta de potencia con igualdad).
  4. ¿La corrección es transferible? Factores Final/Late por celda × mes calendario ajustados en el
     año 1 y evaluados en el año 2 (validación cruzada temporal). Si no transfieren, no sirven para 2026.

    python validate_late.py              # descarga lo que falte (1 job Harmony por mes) y evalúa
    python validate_late.py --evaluate   # sólo evalúa lo ya descargado

Salidas en data/analysis/late_vs_final_*.{csv,json}.
"""
import argparse
import json
import re
import sys
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from harmony import Collection, Request
from scipy import stats

import config
from analysis import OUT_DIR, WET_MM, open_grid, _json_default
from download import login, mexico_geojson, month_chunks
from process import mexico_mask

OVERLAP_DIR = config.BASE_DIR / "data" / "raw_late_overlap"
DATE_RE = re.compile(r"3IMERG\.(\d{8})-")
MARGIN = 0.10          # margen de equivalencia: ±10 % en totales
WET_PP_MARGIN = 5.0    # ±5 pp en frecuencia de días con lluvia


# ---------------------------------------------------------------------------
def download(start: str, end: str) -> None:
    harmony = login()
    shape = mexico_geojson()
    collection = Collection(id=config.RUNS["late"]["concept_id"])
    for t0, t1 in month_chunks(start, end):
        out_dir = OVERLAP_DIR / t0.strftime("%Y-%m")
        expected = (t1.date() - t0.date()).days + 1
        present = {DATE_RE.search(f.name).group(1) for f in out_dir.glob("*.nc*") if DATE_RE.search(f.name)}
        if len(present) >= expected:
            print(f"↷ {t0:%Y-%m} completo, se omite")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        for partial in out_dir.glob("*.nc*"):
            partial.unlink()
        req = Request(collection=collection, shape=str(shape), temporal={"start": t0, "stop": t1},
                      variables=[config.VARIABLE])
        if not req.is_valid():
            sys.exit(f"Petición inválida: {req.error_messages()}")
        job = harmony.submit(req)
        print(f"→ {t0:%Y-%m}: job Harmony {job}")
        harmony.wait_for_processing(job, show_progress=False)
        files = [f.result() for f in harmony.download_all(job, directory=str(out_dir), overwrite=True)]
        print(f"✔ {t0:%Y-%m}: {len(files)} archivos")


def open_late() -> xr.DataArray:
    files = sorted(OVERLAP_DIR.rglob("*.nc*"))
    if not files:
        sys.exit(f"No hay archivos en {OVERLAP_DIR}; corre sin --evaluate.")
    by_day = {}
    for f in files:
        m = DATE_RE.search(f.name)
        if m and (m.group(1) not in by_day or f.stat().st_mtime > by_day[m.group(1)].stat().st_mtime):
            by_day[m.group(1)] = f
    ds = xr.open_mfdataset(sorted(by_day.values()), combine="by_coords", engine="h5netcdf")
    pr = ds[config.VARIABLE].transpose("time", "lat", "lon")
    pr["time"] = pr["time"].dt.floor("D")
    return mexico_mask(pr.sortby("lat").sortby("time").load())


# ---------------------------------------------------------------------------
def tost_log_ratio(late: np.ndarray, final: np.ndarray, margin: float = MARGIN) -> dict:
    """Equivalencia de totales emparejados: log(L/F) dentro de ±log(1±margen). Devuelve el p del TOST
    (equivalente si p < 0.05) y la mediana del cociente con IC95 bootstrap."""
    ok = (late > 0) & (final > 0)
    d = np.log(late[ok] / final[ok])
    if len(d) < 5:
        return {"n": int(len(d))}
    lo, hi = np.log(1 - margin), np.log(1 + margin)
    p_lo = stats.ttest_1samp(d, lo, alternative="greater").pvalue
    p_hi = stats.ttest_1samp(d, hi, alternative="less").pvalue
    rng = np.random.default_rng(0)
    boot = [np.median(rng.choice(d, len(d))) for _ in range(2000)]
    return {"n": int(len(d)), "cociente_mediano": round(float(np.exp(np.median(d))), 3),
            "ic95": [round(float(np.exp(np.percentile(boot, 2.5))), 3), round(float(np.exp(np.percentile(boot, 97.5))), 3)],
            "p_tost": float(f"{max(p_lo, p_hi):.2g}"), "equivalente_±10%": bool(max(p_lo, p_hi) < 0.05)}


def detection(late: np.ndarray, final: np.ndarray) -> dict:
    ok = ~np.isnan(late) & ~np.isnan(final)
    L, F = late[ok] >= WET_MM, final[ok] >= WET_MM
    hits, miss, fa = np.sum(L & F), np.sum(~L & F), np.sum(L & ~F)
    return {"frec_final_pct": round(100 * F.mean(), 1), "frec_late_pct": round(100 * L.mean(), 1),
            "dif_pp": round(100 * (L.mean() - F.mean()), 1),
            "POD": round(hits / (hits + miss), 3), "FAR": round(fa / (hits + fa), 3), "CSI": round(hits / (hits + miss + fa), 3)}


def evaluate() -> None:
    L = open_late()
    F = open_grid()["precipitation"].sel(time=L.time).load()          # Final, mismos días y celdas
    t = L.time.to_index()
    print(f"Traslape Late/Final: {t[0].date()} → {t[-1].date()} ({len(t)} días)")
    valid = (F.notnull().any("time") & L.notnull().any("time")).values

    # Regiones climáticas (clusters.nc) para ver si el sesgo depende del régimen
    cl = xr.open_dataset(OUT_DIR / "clusters.nc")["cluster"].values if (OUT_DIR / "clusters.nc").exists() else None
    w = np.cos(np.deg2rad(F.lat.values))[:, None] * np.ones(len(F.lon))

    def regional(da, mask):
        x = da.values[:, mask & valid]; ww = w[mask & valid]
        ok = ~np.isnan(x)
        return (np.nan_to_num(x) * ww).sum(1) / (ok * ww).sum(1)

    regions = {"México": np.ones_like(valid)}
    if cl is not None:
        for c in range(1, int(np.nanmax(cl)) + 1):
            regions[f"grupo {c}"] = cl == c

    # 1) Totales mensuales por celda: cociente L/F por mes calendario
    Lm, Fm = L.resample(time="MS").sum(min_count=28), F.resample(time="MS").sum(min_count=28)
    rows = []
    for tm in Lm.time.to_index():
        l, f = Lm.sel(time=tm).values[valid], Fm.sel(time=tm).values[valid]
        ok = (f > 5) & ~np.isnan(l)          # celdas con lluvia apreciable (evita cocientes de ~0/~0)
        rows.append({"mes": str(tm.date())[:7], "celdas": int(ok.sum()),
                     "cociente_mediano_celda": round(float(np.median(l[ok] / f[ok])), 3) if ok.any() else None,
                     "r_espacial": round(float(np.corrcoef(l[ok], f[ok])[0, 1]), 3) if ok.sum() > 10 else None,
                     **{f"{k}_L/F": round(float(regional(Lm.sel(time=[tm]), m)[0] / regional(Fm.sel(time=[tm]), m)[0]), 3)
                        for k, m in regions.items()}})
    monthly = pd.DataFrame(rows).set_index("mes")

    # 2) Temporadas (oficiales, corte diario) por región
    from analysis import in_season
    seas_rows, tost_rows = [], {}
    idx = t
    for season in ("lluvias", "secas"):
        mask_t = in_season(idx, season)
        key = idx.year + ((idx.month >= 11).astype(int) if season == "secas" else 0)
        for k_name, m in regions.items():
            lr, fr = regional(L, m), regional(F, m)
            for y in sorted(set(key[mask_t])):
                sel = mask_t & (key == y)
                window = (pd.date_range(f"{y - 1}-11-16", f"{y}-05-14") if season == "secas"
                          else pd.date_range(f"{y}-05-15", f"{y}-11-15"))
                if sel.sum() < 0.9 * len(window):          # temporada incompleta en el traslape
                    continue
                seas_rows.append({"region": k_name, "temporada": f"{season} {y}", "final_mm": round(float(fr[sel].sum()), 1),
                                  "late_mm": round(float(lr[sel].sum()), 1),
                                  "dif_pct": round(float(100 * (lr[sel].sum() / fr[sel].sum() - 1)), 1)})
    seasons = pd.DataFrame(seas_rows)

    # 3) Equivalencia (TOST) de totales REGIONALES por mes (región climática × mes). La unidad es la región,
    #    no la celda: las celdas vecinas comparten el clima del día, así que 200 mil celdas-mes NO son
    #    independientes y un TOST por celda "pasa" con cualquier sesgo medio pequeño.
    reg_cols = [c for c in monthly.columns if c.endswith("_L/F")]
    for season, months in (("lluvias (may–nov)", [5, 6, 7, 8, 9, 10, 11]), ("secas (dic–abr)", [12, 1, 2, 3, 4])):
        sel = pd.to_datetime(monthly.index).month.isin(months)
        r = monthly.loc[sel, reg_cols].values.ravel().astype(float)
        r = r[np.isfinite(r) & (r > 0)]
        tost_rows[season] = tost_log_ratio(r, np.ones_like(r))      # cociente L/F vs 1
        tost_rows[season]["unidad"] = "región × mes"

    # 4) Detección de días con lluvia
    det = {"México": detection(L.values[:, valid].ravel(), F.values[:, valid].ravel())}
    for k_name, m in list(regions.items())[1:]:
        det[k_name] = detection(L.values[:, m & valid].ravel(), F.values[:, m & valid].ravel())

    # 5) ¿La corrección transfiere entre años? Factor F/L por celda × mes calendario (año 1 → año 2)
    yr1 = (Lm.time.to_index() < pd.Timestamp(t[0]) + pd.DateOffset(years=1))
    transfer = {}
    if (~yr1).sum() >= 12:
        fac = (Fm.isel(time=yr1).groupby("time.month").sum() / Lm.isel(time=yr1).groupby("time.month").sum())
        fac = fac.where(np.isfinite(fac)).clip(0.2, 5)
        L2, F2 = Lm.isel(time=~yr1), Fm.isel(time=~yr1)
        corr = L2 * fac.sel(month=L2.time.dt.month).drop_vars("month")
        for name, est in (("sin_corregir", L2), ("corregido_año1", corr)):
            l, f = est.values[:, valid].ravel(), F2.values[:, valid].ravel()
            ok = (f > 5) & ~np.isnan(l)
            r = {"cociente_mediano_celda": round(float(np.median(l[ok] / f[ok])), 3)}   # descriptivo (sin TOST por celda)
            r["error_abs_mediano_pct"] = round(float(100 * np.median(np.abs(l[ok] / f[ok] - 1))), 1)
            nat_l = regional(est, regions["México"]).sum(); nat_f = regional(F2, regions["México"]).sum()
            r["sesgo_nacional_año2_pct"] = round(float(100 * (nat_l / nat_f - 1)), 1)
            transfer[name] = r

    # 6) Contexto: ¿el sesgo es pequeño frente a la variabilidad interanual? (mensual 1998–2024)
    zm = xr.open_dataset(config.ZONAL_MONTHLY)
    nat = zm.precipitation.sel(region="00").to_series()
    nat = nat[nat.index.year <= 2024]
    seasonal = nat[nat.index.month.isin([5, 6, 7, 8, 9, 10, 11])].groupby(lambda x: x.year).sum()
    cv_interanual = float(seasonal.std() / seasonal.mean())

    out = {"traslape": [str(t[0].date()), str(t[-1].date())], "margen_equivalencia": MARGIN,
           "equivalencia_totales_mensuales": tost_rows, "deteccion_dias_lluvia": det,
           "transferencia_correccion": transfer,
           "cv_interanual_temporada_lluvias_nacional": round(cv_interanual, 3)}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    monthly.to_csv(OUT_DIR / "late_vs_final_mensual.csv")
    seasons.to_csv(OUT_DIR / "late_vs_final_temporadas.csv", index=False)
    (OUT_DIR / "late_vs_final_resumen.json").write_text(json.dumps(out, ensure_ascii=False, indent=2, default=_json_default))

    pd.set_option("display.width", 200)
    print("\n== Cociente Late/Final de totales mensuales (regional, ponderado) ==")
    print(monthly.to_string())
    print("\n== Temporadas (corte oficial) ==")
    print(seasons.pivot_table(index="region", columns="temporada", values="dif_pct").to_string()
          if len(seasons) else "  (ninguna temporada completa en el traslape)")
    print("\n== Resumen ==")
    print(json.dumps(out, ensure_ascii=False, indent=2, default=_json_default))
    print("  → data/analysis/late_vs_final_{mensual.csv, temporadas.csv, resumen.json}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2023-10-01")
    ap.add_argument("--end", default="2025-09-30")
    ap.add_argument("--evaluate", action="store_true")
    a = ap.parse_args()
    if not a.evaluate:
        download(a.start, a.end)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        evaluate()


if __name__ == "__main__":
    main()

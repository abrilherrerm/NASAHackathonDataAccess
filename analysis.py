"""Paso 4: análisis estadístico de la precipitación IMERG.

Herramientas para preguntas tipo "¿cuándo conviene viajar a Puerto Vallarta?" o
"¿cuántas temporadas de lluvia en la ZMVM han rebasado un umbral de anomalía?",
con pruebas propias de precipitación y candados contra errores tipo I y tipo II.

    python analysis.py shape    [--cve 14067 ...]       # forma de los datos -> qué métodos aplican
    python analysis.py travel   [--cve 14067] [--trip 7] # mejor época para viajar (Puerto Vallarta)
    python analysis.py season   [--region zmvm]          # temporada de lluvias vs umbrales de anomalía
    python analysis.py trends   [--cve 09]               # Mann-Kendall estacional / Hamed-Rao / Sen
    python analysis.py power                             # tendencia mínima detectable según nº de años
    python analysis.py clusters                          # regiones homogéneas + zonas de alto gradiente
    python analysis.py field                             # tendencia por celda con control FDR

Los resultados se guardan en data/analysis/. Ver ANALISIS.md para el porqué de cada método.
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pymannkendall as mk
import xarray as xr
from scipy import stats
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import cdist

import config

OUT_DIR = config.BASE_DIR / "data" / "analysis"

# ---------------------------------------------------------------------------
# Umbrales y candados (ver ANALISIS.md §2)
# ---------------------------------------------------------------------------
WET_MM = 1.0              # día con lluvia (umbral OMM/ETCCDI; IMERG tiene "llovizna" espuria < 1 mm)
HEAVY_MM = 10.0           # día de lluvia fuerte (ETCCDI R10mm)
VERY_HEAVY_MM = 20.0      # ETCCDI R20mm
ALPHA = 0.05              # nivel de significancia por prueba
FDR_Q = 0.10              # tasa de falsos descubrimientos para campos (Wilks 2016)
MIN_DAY_FRAC = 0.9        # un mes cuenta como completo con >= 90 % de días
MIN_YEARS_TREND = 10      # por debajo, Mann-Kendall no tiene potencia (ver `power`)
MIN_YEARS_BASELINE = 20   # climatología de referencia (la OMM recomienda 30)
MIN_YEARS_YEAR_BOOT = 8   # bootstrap por años (captura variabilidad interanual / ENSO)

# Regiones predefinidas (claves INEGI CVEGEO). Editables.
REGIONS = {
    "pv": ["14067"],                      # Puerto Vallarta
    "pv_bahia": ["14067", "18020"],       # + Bahía de Banderas (Nayarit), misma bahía
    # Núcleo conurbado de la ZMVM: las 16 alcaldías (entidad 09) + 17 municipios del Edomex.
    # Aproximación: la delimitación oficial (SEDATU/CONAPO/INEGI) tiene 76 municipios.
    "zmvm": ["09", "15033", "15058", "15057", "15104", "15031", "15013", "15121", "15109",
             "15039", "15122", "15037", "15020", "15070", "15081", "15060", "15099", "15025"],
}
REGION_LABELS = {"pv": "Puerto Vallarta", "pv_bahia": "Bahía de Banderas (PV + Nayarit)",
                 "zmvm": "ZMVM núcleo (CDMX + 17 mun. Edomex)"}
MONTHS_ES = ["ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"]


def warn(msg: str):
    print(f"  ⚠ {msg}")


def save(name: str, obj):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    if isinstance(obj, pd.DataFrame):
        obj.to_csv(path, index=True)
    else:
        path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=_json_default))
    print(f"  → {path.relative_to(config.BASE_DIR)}")


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (pd.Timestamp,)):
        return str(o.date())
    raise TypeError(type(o))


# ---------------------------------------------------------------------------
# Carga de series
# ---------------------------------------------------------------------------
def load_zonal() -> xr.Dataset:
    if not config.ZONAL_SERIES.exists():
        sys.exit("No existe zonal_series.nc. Corre make_boundaries.py y zonal.py.")
    return xr.open_dataset(config.ZONAL_SERIES).load()


def region_series(cves: list[str], zs: xr.Dataset | None = None) -> tuple[pd.Series, str]:
    """Serie diaria (mm/día) de una o varias regiones INEGI, promediada por área.
    Los días faltantes en el archivo quedan como NaN (calendario diario completo)."""
    zs = load_zonal() if zs is None else zs
    missing = [c for c in cves if c not in zs.region.values]
    if missing:
        sys.exit(f"Claves no encontradas: {missing}")
    sub = zs.sel(region=cves)
    w = sub.area_km2
    da = (sub.precipitation * w).sum("region", min_count=1) / w.sum()
    s = da.to_series()
    s = s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="D"))
    names = [str(n) for n in sub.nombre.values]
    label = names[0] if len(names) == 1 else f"{names[0]} + {len(names) - 1} regiones"
    return s.astype("float64"), label


def point_series(lat: float, lon: float) -> tuple[pd.Series, str]:
    pr = open_grid()["precipitation"]
    # Celda con dato más cercana: en la costa la celda exacta puede estar enmascarada como mar
    valid = pr.isel(time=slice(0, 31)).notnull().any("time").values
    LA, LO = np.meshgrid(pr.lat.values, pr.lon.values, indexing="ij")
    d2 = np.where(valid, (LA - lat) ** 2 + ((LO - lon) * np.cos(np.deg2rad(lat))) ** 2, np.inf)
    i, j = np.unravel_index(np.argmin(d2), d2.shape)
    da = pr.isel(lat=i, lon=j).load()
    s = da.to_series()
    s = s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="D"))
    return s.astype("float64"), f"Celda ({float(da.lat):.2f}, {float(da.lon):.2f})"


def open_grid() -> xr.Dataset:
    files = sorted(config.PROCESSED_DIR.glob("imerg_mexico_daily_*.nc"))
    if not files:
        sys.exit("No hay NetCDF procesados. Corre process.py.")
    return xr.open_mfdataset(files, combine="by_coords", chunks={})


def monthly_totals(s: pd.Series) -> pd.Series:
    """Totales mensuales; NaN si el mes tiene < MIN_DAY_FRAC de días con dato."""
    tot = s.resample("MS").sum(min_count=1)
    frac = s.notna().resample("MS").sum() / s.index.to_series().resample("MS").count()
    return tot.where(frac >= MIN_DAY_FRAC)


def complete_years(monthly: pd.Series, months=range(1, 13)) -> list[int]:
    """Años en que todos los `months` están completos."""
    months = list(months)
    m = monthly[monthly.index.month.isin(months)]
    ok = m.groupby(m.index.year).apply(lambda g: g.notna().sum() == len(months))
    return [int(y) for y in ok[ok].index]


def resolve_series(args) -> tuple[pd.Series, str]:
    if getattr(args, "point", None):
        return point_series(*args.point)
    if getattr(args, "region", None):
        s, _ = region_series(REGIONS[args.region])
        return s, REGION_LABELS[args.region]
    return region_series(args.cve)


# ---------------------------------------------------------------------------
# 1. Forma de los datos -> qué métodos aplican
# ---------------------------------------------------------------------------
def best_distribution(x: np.ndarray) -> dict:
    """Compara gamma / lognormal / Weibull (loc=0) por AIC sobre valores positivos."""
    x = x[x > 0]
    if len(x) < 30:
        return {"n": len(x), "best": None}
    out = {}
    for name, dist in [("gamma", stats.gamma), ("lognorm", stats.lognorm), ("weibull", stats.weibull_min)]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            params = dist.fit(x, floc=0)
        ll = np.sum(dist.logpdf(x, *params))
        out[name] = round(2 * (len(params) - 1) - 2 * ll, 1)   # loc fijo: no cuenta como parámetro
    return {"n": len(x), "aic": out, "best": min(out, key=out.get)}


def shape_report(s: pd.Series, label: str) -> dict:
    d = s.dropna()
    m = monthly_totals(s).dropna()
    clim = m.groupby(m.index.month).transform("mean")
    anom = m - clim
    n_per_month = m.groupby(m.index.month).size()
    years = complete_years(monthly_totals(s))
    lag1_d = d.autocorr(1)
    lag1_anom = anom.autocorr(1) if len(anom) > 12 else np.nan
    # Significancia aproximada de r1 (Bartlett: ±1.96/√n)
    lag1_anom_sig = bool(abs(lag1_anom) > 1.96 / np.sqrt(len(anom))) if len(anom) > 12 else None
    # Estacionalidad: Kruskal-Wallis de totales mensuales por mes calendario
    groups = [g.values for _, g in m.groupby(m.index.month) if len(g) >= 2]
    kw_p = stats.kruskal(*groups).pvalue if len(groups) >= 3 else np.nan

    rep = {
        "serie": label,
        "periodo": f"{d.index.min().date()} → {d.index.max().date()}",
        "dias_con_dato": int(d.size),
        "anios_completos": years,
        "muestras_por_mes_calendario": {MONTHS_ES[k - 1]: int(v) for k, v in n_per_month.items()},
        "frac_dias_secos(<0.1mm)": round(float((d < 0.1).mean()), 3),
        "frac_dias_lluvia(>=1mm)": round(float((d >= WET_MM).mean()), 3),
        "asimetria_diaria": round(float(stats.skew(d)), 2),
        "autocorr_lag1_diaria": round(float(lag1_d), 2),
        "autocorr_lag1_anomalias_mensuales": None if np.isnan(lag1_anom) else round(float(lag1_anom), 2),
        "estacionalidad_kruskal_p": None if np.isnan(kw_p) else float(f"{kw_p:.2g}"),
        "distribucion_dias_lluvia": best_distribution(d[d >= WET_MM].values),
        "distribucion_totales_mensuales": best_distribution(m.values),
    }
    # Recomendaciones derivadas de la forma de los datos
    rec = []
    if rep["frac_dias_secos(<0.1mm)"] > 0.3 or abs(rep["asimetria_diaria"]) > 1:
        rec.append("Inflada en ceros y asimétrica: usar pruebas no paramétricas (Mann-Kendall, Sen, "
                   "Kruskal-Wallis, bootstrap) o modelos gamma con masa en cero (SPI). Evitar t-test/OLS sobre mm.")
    if not np.isnan(kw_p) and kw_p < ALPHA:
        rec.append("Estacionalidad significativa: comparar contra la climatología del mismo mes y usar "
                   "Mann-Kendall estacional (Hirsch 1982) para tendencias mensuales.")
    if lag1_d > 0.1:
        rec.append(f"Persistencia diaria (r1={lag1_d:.2f}): los días no son independientes -> bootstrap por "
                   "bloques para IC de frecuencias diarias.")
    if lag1_anom_sig:
        rec.append("Autocorrelación en anomalías mensuales: usar Hamed-Rao (o Mann-Kendall estacional "
                   "correlacionado) para no inflar el error tipo I.")
    if len(years) < MIN_YEARS_TREND:
        rec.append(f"Sólo {len(years)} años completos: NO hay potencia para tendencias "
                   f"(mínimo {MIN_YEARS_TREND}; ver `python analysis.py power`). Ampliar la descarga "
                   "(IMERG V07 diario existe desde 1998).")
    if len(years) < MIN_YEARS_BASELINE:
        rec.append(f"Sin climatología de referencia confiable (< {MIN_YEARS_BASELINE} años): las anomalías y "
                   "umbrales (SPI, percentiles) serían autorreferidos.")
    dist = rep["distribucion_dias_lluvia"]["best"]
    if dist:
        rec.append(f"Montos de días con lluvia: mejor ajuste = {dist} (AIC).")
    rep["recomendaciones"] = rec
    return rep


def cmd_shape(args):
    s, label = resolve_series(args)
    rep = shape_report(s, label)
    print(json.dumps(rep, ensure_ascii=False, indent=2, default=_json_default))
    save(f"shape_{_slug(label)}.json", rep)


# ---------------------------------------------------------------------------
# 2. Bootstrap respetando la dependencia temporal
# ---------------------------------------------------------------------------
def _block_resample(segments: list[np.ndarray], block: int, rng) -> np.ndarray:
    """Moving block bootstrap dentro de segmentos contiguos (no cruza de un año a otro)."""
    pool = [seg[i:i + block] for seg in segments for i in range(max(1, len(seg) - block + 1))]
    total = sum(len(seg) for seg in segments)
    out, n = [], 0
    while n < total:
        b = pool[rng.integers(len(pool))]
        out.append(b); n += len(b)
    return np.concatenate(out)[:total]


def bootstrap_stat(segments: list[np.ndarray], fn, block: int, n_boot: int, rng, by_year: bool) -> np.ndarray:
    """Distribución bootstrap de fn(valores).
    by_year=True: remuestrea años completos (captura variabilidad interanual; requiere muchos años).
    by_year=False: bloques dentro de cada año (sólo variabilidad intraestacional -> IC optimistas)."""
    res = np.empty(n_boot)
    for b in range(n_boot):
        if by_year:
            pick = rng.integers(len(segments), size=len(segments))
            res[b] = fn(np.concatenate([segments[i] for i in pick]))
        else:
            res[b] = fn(_block_resample(segments, block, rng))
    return res


def holm(pvals: dict) -> dict:
    """Corrección de Holm-Bonferroni (controla FWER en comparaciones múltiples)."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, running = len(items), {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        out[k] = running
    return out


def bh_fdr(p: np.ndarray, q: float) -> np.ndarray:
    """Benjamini-Hochberg: máscara de descubrimientos con FDR <= q (ignora NaN)."""
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    pv = p[ok]
    order = np.argsort(pv)
    m = len(pv)
    thresh = q * (np.arange(1, m + 1) / m)
    passed = pv[order] <= thresh
    k = np.max(np.nonzero(passed)[0]) + 1 if passed.any() else 0
    sig = np.zeros(m, bool); sig[order[:k]] = True
    out = np.zeros(p.shape, bool); out[ok] = sig
    return out


# ---------------------------------------------------------------------------
# 3. Mejor época para viajar
# ---------------------------------------------------------------------------
def cmd_travel(args):
    s, label = resolve_series(args)
    trip, rng = args.trip, np.random.default_rng(args.seed)
    wet = (s >= WET_MM).astype(float).where(s.notna())
    heavy = (s >= HEAVY_MM).astype(float).where(s.notna())
    # "Viaje sin lluvia": ventana de `trip` días que empieza ese día con <= 1 día de lluvia
    wet_in_trip = wet[::-1].rolling(trip, min_periods=trip).sum()[::-1]
    good = (wet_in_trip <= 1).astype(float).where(wet_in_trip.notna())

    years = sorted(set(s.dropna().index.year))
    n_years_by_month = s.dropna().groupby(s.dropna().index.month).apply(lambda g: g.index.year.nunique())
    by_year = int(n_years_by_month.min()) >= MIN_YEARS_YEAR_BOOT
    print(f"\n{label} · viaje de {trip} días · años con dato: {years[0]}–{years[-1]}")
    if not by_year:
        warn(f"Sólo {int(n_years_by_month.min())}–{int(n_years_by_month.max())} años por mes: IC por bloques "
             "dentro de cada año. NO incluyen variabilidad interanual (ENSO, huracanes): son optimistas. "
             f"Con >= {MIN_YEARS_YEAR_BOOT} años se usa bootstrap por años automáticamente.")

    def segs(x: pd.Series, month: int) -> list[np.ndarray]:
        x = x[x.index.month == month].dropna()
        return [g.values for _, g in x.groupby(x.index.year) if len(g)]

    rows, boots = [], {}
    for mth in range(1, 13):
        sw, sg, sh = segs(wet, mth), segs(good, mth), segs(heavy, mth)
        if not sw:
            continue
        bw = bootstrap_stat(sw, np.mean, 7, args.boot, rng, by_year)
        bg = bootstrap_stat(sg, np.mean, 2 * trip, args.boot, rng, by_year) if sg else np.full(args.boot, np.nan)
        boots[mth] = bg if args.rank_by == "good" else bw
        mm = s[s.index.month == mth].dropna()
        rows.append({
            "mes": MONTHS_ES[mth - 1], "n_dias": len(mm), "n_anios": mm.index.year.nunique(),
            "mm_dia_medio": round(mm.mean(), 2),
            "p_lluvia": round(np.concatenate(sw).mean(), 3),
            "p_lluvia_ic95": f"{np.nanpercentile(bw, 2.5):.3f}–{np.nanpercentile(bw, 97.5):.3f}",
            "p_lluvia_fuerte": round(np.concatenate(sh).mean(), 3),
            "dias_lluvia_esperados_viaje": round(trip * np.concatenate(sw).mean(), 1),
            "p_viaje_seco": round(np.concatenate(sg).mean(), 3) if sg else np.nan,
            "p_viaje_seco_ic95": f"{np.nanpercentile(bg, 2.5):.3f}–{np.nanpercentile(bg, 97.5):.3f}" if sg else "",
        })
    df = pd.DataFrame(rows).set_index("mes")
    key = "p_viaje_seco" if args.rank_by == "good" else "p_lluvia"
    best_m = df[key].idxmax() if args.rank_by == "good" else df[key].idxmin()
    best_i = MONTHS_ES.index(best_m) + 1

    # ¿Qué meses son estadísticamente indistinguibles del mejor? (bootstrap de diferencias + Holm)
    raw_p, h = {}, {}
    for mth, b in boots.items():
        if mth == best_i:
            continue
        diff = b - boots[best_i]
        raw_p[mth] = min(1.0, 2 * min(np.mean(diff <= 0), np.mean(diff >= 0)) + 1 / args.boot)
        p1, p2 = df.loc[MONTHS_ES[mth - 1], key], df.loc[best_m, key]
        h[mth] = round(abs(2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p2))), 2)   # h de Cohen
    adj = holm(raw_p)
    df["p_vs_mejor_holm"] = [np.nan if (i := MONTHS_ES.index(m) + 1) == best_i else round(adj[i], 4) for m in df.index]
    df["h_cohen_vs_mejor"] = [np.nan if (i := MONTHS_ES.index(m) + 1) == best_i else h[i] for m in df.index]
    df["recomendable"] = [(m == best_m) or (adj[MONTHS_ES.index(m) + 1] >= ALPHA) for m in df.index]

    print(df.to_string())
    window = ", ".join(df.index[df["recomendable"]])
    print(f"\n  Mejor mes por {key}: {best_m}. Indistinguibles del mejor (Holm, α={ALPHA}): {window}")
    print("  h de Cohen: 0.2 pequeño · 0.5 mediano · 0.8 grande (relevancia práctica, no sólo significancia).")
    save(f"travel_{_slug(label)}.csv", df)


# ---------------------------------------------------------------------------
# 4. Temporada de lluvias vs umbrales de anomalía
# ---------------------------------------------------------------------------
def spi(values: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    """SPI (McKee 1993; OMM 2012): gamma ajustada a la referencia con masa en cero."""
    q = np.mean(baseline <= 0)
    pos = baseline[baseline > 0]
    a, _, scale = stats.gamma.fit(pos, floc=0)
    cdf = q + (1 - q) * stats.gamma.cdf(np.clip(values, 0, None), a, scale=scale)
    return stats.norm.ppf(np.clip(cdf, 1e-6, 1 - 1e-6))


def spi_class(x: np.ndarray) -> np.ndarray:
    """Categorías de McKee et al. (1993): <= -1 seco, >= +1 húmedo (mismos cortes que los rebases)."""
    conds = [x <= -2, x <= -1.5, x <= -1, x < 1, x < 1.5, x < 2]
    names = ["extremo seco", "severo seco", "moderado seco", "normal", "moderado húmedo", "severo húmedo"]
    return np.select(conds, names, default="extremo húmedo")


def in_season(index: pd.DatetimeIndex, season: str) -> np.ndarray:
    """Temporadas oficiales: lluvias 15 may–15 nov; secas = complemento exacto (16 nov–14 may)."""
    m, d = index.month, index.day
    wet = ((m == 5) & (d >= 15)) | m.isin([6, 7, 8, 9, 10]) | ((m == 11) & (d <= 15))
    return np.asarray(wet if season == "lluvias" else ~wet)


SEASON_NAMES = {"lluvias": "lluvias 15may–15nov", "secas": "secas 16nov–14may"}


def window_days(s: pd.Series, months=None, season=None) -> tuple[pd.Series, pd.Series, list[int]]:
    """Días dentro de la ventana (temporada oficial o lista de meses), su temporada-año (las que
    cruzan el año cuentan en el año en que terminan) y las temporadas con >= MIN_DAY_FRAC de días."""
    idx = s.index
    if season:
        mask, start = in_season(idx, season), (11 if season == "secas" else None)
    else:
        months = list(months)
        mask, start = np.asarray(idx.month.isin(months)), (months[0] if months != sorted(months) else None)
    key = pd.Series(idx.year + ((idx.month >= start).astype(int) if start else 0), index=idx)
    x, k = s[mask], key[mask]
    frac = x.notna().groupby(k).mean()
    return x, k, [int(y) for y, f in frac.items() if f >= MIN_DAY_FRAC]


def season_indices(s: pd.Series, months, season, p95_wet: float | None) -> pd.DataFrame:
    """Índices ETCCDI por temporada (año) con temporada completa."""
    xs, ks, years = window_days(s, months, season)
    rows = {}
    for y in years:
        x = xs[ks == y]
        dry_runs = (x < WET_MM).astype(int)
        cdd = dry_runs.groupby((dry_runs == 0).cumsum()).sum().max()
        rows[y] = {
            "total_mm": round(x.sum(), 1),
            "dias_lluvia": int((x >= WET_MM).sum()),
            "SDII_mm": round(x[x >= WET_MM].mean(), 2),           # intensidad media en días con lluvia
            "Rx1day_mm": round(x.max(), 1),
            "Rx5day_mm": round(x.rolling(5).sum().max(), 1),
            "R20mm_dias": int((x >= VERY_HEAVY_MM).sum()),
            "R95p_dias": int((x > p95_wet).sum()) if p95_wet else None,
            "CDD_dias": int(cdd),                                  # racha seca más larga en la temporada
        }
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("anio")


def cmd_season(args):
    s, label = resolve_series(args)
    months, season = args.months, (None if args.months else args.window)
    base0, base1 = args.base
    xs, ks, years = window_days(s, months, season)
    base_years = [y for y in years if base0 <= y <= base1]
    seas = SEASON_NAMES[season].split()[0] if season else "".join(MONTHS_ES[m - 1][0].upper() for m in months)
    print(f"\n{label} · temporada {SEASON_NAMES.get(season, seas)} · temporadas completas: {years or 'ninguna'}")

    in_base = xs[ks.isin(base_years)]
    p95 = float(np.nanpercentile(in_base[in_base >= WET_MM], 95)) if len(in_base) else None
    df = season_indices(s, months, season, p95)
    if df.empty:
        warn("No hay ninguna temporada completa."); return

    enough = len(base_years) >= MIN_YEARS_BASELINE
    if not enough:
        warn(f"Referencia {base0}–{base1}: {len(base_years)} temporadas (< {MIN_YEARS_BASELINE}). "
             "Las anomalías serían autorreferidas y el umbral inestable.")
    if not base_years and args.force:
        warn(f"No hay temporadas completas en {base0}–{base1}; ajusta --base a los años disponibles: {years}.")
    if enough or (args.force and base_years):
        tag = "" if enough else " (PROVISIONAL)"
        ref = df.loc[base_years]
        df["anomalia_pct"] = ((df.total_mm / ref.total_mm.median() - 1) * 100).round(1)
        df["SPI"] = spi(df.total_mm.values, ref.total_mm.values).round(2)
        df["Rx1day_>P90ref"] = df.Rx1day_mm > np.percentile(ref.Rx1day_mm, 90)
        df["clase"] = spi_class(df.SPI.values)
        # Frecuencia de rebase en el periodo de evaluación vs lo esperado bajo la referencia
        ev = df.loc[[y for y in df.index if args.eval[0] <= y <= args.eval[1]]]
        summary = {}
        for name, mask, p0 in [("SPI<=-1 (seco)", ev.SPI <= -1, stats.norm.cdf(-1)),
                               ("SPI>=+1 (húmedo)", ev.SPI >= 1, stats.norm.cdf(-1)),
                               ("Rx1day > P90 ref", ev["Rx1day_>P90ref"], 0.10)]:
            k, n = int(mask.sum()), len(ev)
            if n == 0:
                continue
            ci = stats.binomtest(k, n).proportion_ci(method="exact")
            summary[name] = {"k": k, "n": n, "frecuencia": round(k / n, 3),
                             "ic95": [round(ci.low, 3), round(ci.high, 3)], "esperado_ref": round(p0, 3),
                             "p_binomial(mayor_que_ref)": round(stats.binomtest(k, n, p0, "greater").pvalue, 4)}
        print(df.to_string())
        print(f"\n  Rebases {args.eval[0]}–{args.eval[1]}{tag}:")
        for k, v in summary.items():
            print(f"    {k}: {v['k']}/{v['n']} = {v['frecuencia']:.0%} (IC95 {v['ic95'][0]:.0%}–{v['ic95'][1]:.0%}; "
                  f"esperado {v['esperado_ref']:.0%}; p={v['p_binomial(mayor_que_ref)']})")
        if ev.shape[0] and ev.shape[0] < 10:
            warn(f"n={ev.shape[0]} temporadas evaluadas: la prueba binomial sólo detectaría cambios enormes (tipo II).")
        save(f"season_{_slug(label)}_{seas}.json", {"serie": label, "temporada": SEASON_NAMES.get(season) or months, "referencia": base_years,
                                                    "provisional": not enough, "p95_dias_lluvia_ref": p95,
                                                    "rebases": summary})
    else:
        print(df.to_string())
        print("  (usa --force para calcular anomalías PROVISIONALES con la referencia disponible)")
    save(f"season_{_slug(label)}_{seas}.csv", df)


# ---------------------------------------------------------------------------
# 5. Tendencias: Mann-Kendall estacional, Hamed-Rao, Sen
# ---------------------------------------------------------------------------
def _mk(res) -> dict:
    return {"tendencia": res.trend, "p": float(f"{res.p:.3g}"), "tau": round(res.Tau, 3),
            "pendiente_sen": round(res.slope, 3)}


def trend_report(s: pd.Series, season_months: list[int]) -> dict:
    mt = monthly_totals(s)
    years = complete_years(mt)
    out = {"anios_completos": years, "n": len(years)}
    if len(years) < MIN_YEARS_TREND:
        out["estado"] = (f"insuficiente: {len(years)} años completos (< {MIN_YEARS_TREND}). Cualquier "
                         "'tendencia' sería ruido o perdería significancia por falta de potencia.")
        return out
    monthly = mt[mt.index.year.isin(years)].values                 # empieza en enero, años completos
    annual = mt[mt.index.year.isin(years)].groupby(lambda t: t.year).sum().values
    seas_years = complete_years(mt, season_months)
    seas = np.array([mt[(mt.index.year == y) & mt.index.month.isin(season_months)].sum() for y in seas_years])

    r1 = pd.Series(annual).autocorr(1)
    out.update({
        "estado": "ok",
        "mensual_MK_estacional": _mk(mk.seasonal_test(monthly, period=12, alpha=ALPHA)),
        "mensual_MK_estacional_correlacionado": _mk(mk.correlated_seasonal_test(monthly, period=12, alpha=ALPHA)),
        "anual_MK_original": _mk(mk.original_test(annual, alpha=ALPHA)),
        "anual_Hamed_Rao": _mk(mk.hamed_rao_modification_test(annual, alpha=ALPHA)),
        "temporada_Hamed_Rao": _mk(mk.hamed_rao_modification_test(seas, alpha=ALPHA)),
        "anual_autocorr_lag1": round(float(r1), 2),
        "unidades_pendiente": "mm/año (anual y temporada) · mm/mes por año (estacional)",
    })
    a, h = out["anual_MK_original"], out["anual_Hamed_Rao"]
    if (a["p"] < ALPHA) != (h["p"] < ALPHA):
        out["nota"] = "MK original y Hamed-Rao discrepan: la autocorrelación cambia la conclusión; reportar Hamed-Rao."
    return out


def cmd_trends(args):
    s, label = resolve_series(args)
    rep = {"serie": label, **trend_report(s, args.months)}
    print(json.dumps(rep, ensure_ascii=False, indent=2, default=_json_default))
    save(f"trends_{_slug(label)}.json", rep)


# ---------------------------------------------------------------------------
# 6. Potencia: tendencia mínima detectable (error tipo II)
# ---------------------------------------------------------------------------
def simulate_power(n: int, beta: float, cv: float, phi: float, sims: int, rng) -> float:
    """Potencia de Hamed-Rao para una tendencia relativa `beta` (fracción de la media por año)
    en totales anuales gamma (coef. de variación `cv`) con AR(1) `phi` (cópula gaussiana)."""
    k = 1 / cv ** 2
    t = np.arange(n) - (n - 1) / 2
    hits = 0
    for _ in range(sims):
        z = np.empty(n); z[0] = rng.standard_normal()
        for i in range(1, n):
            z[i] = phi * z[i - 1] + np.sqrt(1 - phi ** 2) * rng.standard_normal()
        x = stats.gamma.ppf(stats.norm.cdf(z), k, scale=1 / k) * np.clip(1 + beta * t, 0.05, None)
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")   # con n muy chico la varianza corregida sale negativa -> p=NaN
            r = mk.hamed_rao_modification_test(x, alpha=ALPHA)
        hits += bool(r.p < ALPHA) and (r.trend == "increasing")
    return hits / sims


def mde(n: int, cv: float, phi: float, sims: int, rng, power: float = 0.8) -> float:
    lo, hi = 0.0, 0.25
    if simulate_power(n, hi, cv, phi, sims, rng) < power:
        return np.inf
    for _ in range(7):
        mid = (lo + hi) / 2
        if simulate_power(n, mid, cv, phi, sims, rng) >= power:
            hi = mid
        else:
            lo = mid
    return hi


def cmd_power(args):
    rng = np.random.default_rng(args.seed)
    print(f"\nTendencia mínima detectable (potencia 80 %, α={ALPHA}, Hamed-Rao, totales anuales gamma)")
    print("  expresada en % de la media por década. CV típico de totales anuales en México: 0.15–0.35.")
    rows = []
    for cv in args.cv:
        for phi in args.phi:
            row = {"CV": cv, "AR1": phi}
            for n in args.n:
                b = mde(n, cv, phi, args.sims, rng)
                row[f"n={n}"] = "no detectable" if np.isinf(b) else f"{b * 1000:.0f} %/déc"
            rows.append(row)
            print("  " + "  ".join(f"{k}: {v}" for k, v in row.items()))
    save("power_mde.csv", pd.DataFrame(rows).set_index(["CV", "AR1"]))


# ---------------------------------------------------------------------------
# 7. Regionalización: grupos con régimen homogéneo + zonas de alto gradiente
# ---------------------------------------------------------------------------
def grid_features() -> tuple[xr.Dataset, np.ndarray, np.ndarray]:
    """Rasgos por celda a partir de la climatología mensual (forma del ciclo anual + magnitud)."""
    pr = open_grid()["precipitation"].load()
    pr = pr.reindex(time=pd.date_range(pr.time.values[0], pr.time.values[-1], freq="D"))
    ndays = pr.notnull().resample(time="MS").sum()
    mon = pr.resample(time="MS").sum(min_count=1)
    expected = xr.DataArray(mon.time.to_index().days_in_month.values, coords={"time": mon.time}, dims="time")
    mon = mon.where(ndays.max(("lat", "lon")) >= MIN_DAY_FRAC * expected)   # sólo meses completos
    clim = mon.groupby("time.month").mean()                                 # (12, lat, lon)
    annual = clim.sum("month", min_count=12)
    wetfreq = (pr >= WET_MM).where(pr.notnull()).mean("time")
    sdii = pr.where(pr >= WET_MM).mean("time")
    shape = clim / annual                                                   # fracción del total por mes
    si = (abs(shape - 1 / 12)).sum("month")                                 # índice de estacionalidad (Walsh & Lawler)
    feats = xr.Dataset({"annual": annual, "wetfreq": wetfreq, "sdii": sdii, "si": si, "shape": shape})
    valid = (annual.notnull() & (annual > 0)).values
    X = np.column_stack([
        *(shape.values[m][valid] for m in range(12)),
        np.log1p(annual.values[valid]), wetfreq.values[valid], sdii.values[valid], si.values[valid],
    ])
    X = (X - X.mean(0)) / X.std(0)
    X[:, :12] *= np.sqrt(4 / 12)   # que los 12 meses pesen como 4 rasgos, no como 12
    return feats, valid, X


def silhouette(X: np.ndarray, labels: np.ndarray) -> float:
    D = cdist(X, X)
    s = []
    for i in range(len(X)):
        same = labels == labels[i]
        if same.sum() <= 1:
            s.append(0.0); continue
        a = D[i, same].sum() / (same.sum() - 1)
        b = min(D[i, labels == k].mean() for k in np.unique(labels) if k != labels[i])
        s.append((b - a) / max(a, b))
    return float(np.mean(s))


def cmd_clusters(args):
    rng = np.random.default_rng(args.seed)
    feats, valid, X = grid_features()
    sample = rng.choice(len(X), size=min(3000, len(X)), replace=False)
    scores, fits = {}, {}
    for k in range(args.kmin, args.kmax + 1):
        cent, lab = kmeans2(X, k, minit="++", seed=args.seed, iter=50)
        fits[k] = lab
        scores[k] = round(silhouette(X[sample], lab[sample]), 3)
        print(f"  k={k}: silueta={scores[k]}")
    k = args.k or max(scores, key=scores.get)
    lab = fits[k]
    # Reetiquetar por precipitación anual media (1 = más seco)
    ann = feats.annual.values[valid]
    order = np.argsort([ann[lab == c].mean() for c in range(k)])
    lab = np.argsort(order)[lab] + 1

    grid = np.full(valid.shape, np.nan); grid[valid] = lab
    # Alto gradiente: decil superior de |∇ log(1 + total anual)| (montaña-costa, p. ej. Bahía de Banderas)
    la = np.log1p(feats.annual.values)
    gy, gx = np.gradient(la, feats.lat.values, feats.lon.values)
    gmag = np.hypot(gy, gx * np.cos(np.deg2rad(feats.lat.values))[:, None] ** -1)
    thr = np.nanpercentile(gmag[valid], 90)
    high_grad = (gmag >= thr) & valid
    # Celdas frontera: algún vecino en otro grupo
    pad = np.pad(grid, 1, constant_values=np.nan)
    neigh = np.stack([pad[1:-1, :-2], pad[1:-1, 2:], pad[:-2, 1:-1], pad[2:, 1:-1]])
    border = valid & np.any((neigh != grid) & ~np.isnan(neigh), axis=0)

    out = xr.Dataset({
        "cluster": (("lat", "lon"), grid), "alto_gradiente": (("lat", "lon"), high_grad),
        "frontera": (("lat", "lon"), border), "grad_log_anual": (("lat", "lon"), gmag),
        "anual_mm": feats.annual, "frac_dias_lluvia": feats.wetfreq, "SDII": feats.sdii, "estacionalidad": feats.si,
    }, coords={"lat": feats.lat, "lon": feats.lon},
       attrs={"k": k, "silueta": json.dumps(scores), "nota": "Rasgos de climatología: requieren años completos."})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(OUT_DIR / "clusters.nc")
    print(f"  → data/analysis/clusters.nc (k={k})")

    rows = []
    for c in range(1, k + 1):
        m = grid == c
        cyc = feats.shape.values[:, m].mean(1)
        rows.append({"grupo": c, "celdas": int(m.sum()), "anual_mm": round(float(np.nanmean(feats.annual.values[m])), 0),
                     "frac_dias_lluvia": round(float(np.nanmean(feats.wetfreq.values[m])), 3),
                     "mes_pico": MONTHS_ES[int(np.argmax(cyc))],
                     "frac_JJAS": round(float(cyc[5:9].sum()), 2),
                     "estacionalidad": round(float(np.nanmean(feats.si.values[m])), 2),
                     "pct_alto_gradiente": round(100 * float(high_grad[m].mean()), 1),
                     "pct_frontera": round(100 * float(border[m].mean()), 1)})
    df = pd.DataFrame(rows).set_index("grupo")
    print(df.to_string())
    for name, (lat, lon) in {"Puerto Vallarta": (20.62, -105.23), "CDMX": (19.43, -99.13)}.items():
        # celda válida más cercana (la costa puede caer en una celda enmascarada como mar)
        LA, LO = np.meshgrid(feats.lat.values, feats.lon.values, indexing="ij")
        d2 = np.where(valid, (LA - lat) ** 2 + ((LO - lon) * np.cos(np.deg2rad(lat))) ** 2, np.inf)
        i, j = np.unravel_index(np.argmin(d2), d2.shape)
        print(f"  {name} ({LA[i, j]:.2f}, {LO[i, j]:.2f}): grupo {grid[i, j]:.0f} · alto gradiente={bool(high_grad[i, j])} · frontera={bool(border[i, j])}")
    print("  En celdas de alto gradiente/frontera, la media del grupo no representa al sitio: analizar la celda o municipio.")
    save("clusters.csv", df)


# ---------------------------------------------------------------------------
# 8. Campo de tendencias con control de FDR (multiplicidad espacial)
# ---------------------------------------------------------------------------
def cmd_field(args):
    pr = open_grid()["precipitation"]
    full = pd.date_range(pd.Timestamp(pr.time.values[0]), pd.Timestamp(pr.time.values[-1]), freq="D")
    days = pd.Series(1, index=pd.DatetimeIndex(pr.time.values)).reindex(full).fillna(0)
    yr_frac = days.groupby(days.index.year).mean()
    years = [int(y) for y, f in yr_frac.items() if f >= MIN_DAY_FRAC and (days.index.year == y).sum() >= 365]
    if len(years) < MIN_YEARS_TREND:
        warn(f"{len(years)} años completos (< {MIN_YEARS_TREND}): campo de tendencias no calculado."); return
    ann = pr.sel(time=pr.time.dt.year.isin(years)).groupby("time.year").sum(min_count=300).load()
    p = np.full(ann.shape[1:], np.nan); slope = p.copy()
    for i in range(ann.shape[1]):
        for j in range(ann.shape[2]):
            x = ann.values[:, i, j]
            if np.isnan(x).any():
                continue
            r = mk.hamed_rao_modification_test(x, alpha=ALPHA)
            p[i, j], slope[i, j] = r.p, r.slope
    sig_raw = (p < ALPHA)
    sig_fdr = bh_fdr(p.ravel(), FDR_Q).reshape(p.shape)
    n = np.sum(~np.isnan(p))
    print(f"  Celdas: {n} · p<{ALPHA} sin corregir: {sig_raw.sum()} ({sig_raw.sum() / n:.1%}, "
          f"~{ALPHA:.0%} esperado por azar) · significativas con FDR q={FDR_Q}: {sig_fdr.sum()}")
    xr.Dataset({"p_hamed_rao": (("lat", "lon"), p), "sen_mm_anio": (("lat", "lon"), slope),
                "significativo_fdr": (("lat", "lon"), sig_fdr)},
               coords={"lat": ann.lat, "lon": ann.lon}).to_netcdf(OUT_DIR / "field_trends.nc")
    print("  → data/analysis/field_trends.nc")


# ---------------------------------------------------------------------------
def _slug(s: str) -> str:
    import unicodedata
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()
    return "".join(c if c.isalnum() else "_" for c in s.lower()).strip("_")[:40]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def target(p, default_cve):
        p.add_argument("--cve", nargs="+", default=default_cve, help="claves INEGI (estado, municipio o '00')")
        p.add_argument("--region", choices=REGIONS, help="región predefinida (sustituye a --cve)")
        p.add_argument("--point", nargs=2, type=float, metavar=("LAT", "LON"), help="celda IMERG más cercana")

    p = sub.add_parser("shape"); target(p, ["14067"]); p.set_defaults(fn=cmd_shape)
    p = sub.add_parser("travel"); target(p, ["14067"])
    p.add_argument("--trip", type=int, default=7, help="duración del viaje en días")
    p.add_argument("--rank-by", choices=["good", "wet"], default="good",
                   help="good = prob. de viaje con <= 1 día de lluvia; wet = prob. diaria de lluvia")
    p.add_argument("--boot", type=int, default=2000); p.add_argument("--seed", type=int, default=0)
    p.set_defaults(fn=cmd_travel)
    p = sub.add_parser("season"); target(p, ["09"])
    p.add_argument("--window", choices=["lluvias", "secas"], default="lluvias",
                   help="temporada oficial: lluvias 15 may–15 nov (def.) o secas 16 nov–14 may")
    p.add_argument("--months", nargs="+", type=int, help="meses personalizados (sustituye a --window), p. ej. 6 7 8 9")
    p.add_argument("--base", nargs=2, type=int, default=[2001, 2020], help="años de referencia")
    p.add_argument("--eval", nargs=2, type=int, default=[2021, 2025], help="años a evaluar contra la referencia")
    p.add_argument("--force", action="store_true", help="anomalías provisionales aunque falten años")
    p.set_defaults(fn=cmd_season)
    p = sub.add_parser("trends"); target(p, ["09"])
    p.add_argument("--months", nargs="+", type=int, default=[5, 6, 7, 8, 9, 10, 11],
                   help="meses de la temporada (def. may–nov: aproximación mensual de la temporada oficial)")
    p.set_defaults(fn=cmd_trends)
    p = sub.add_parser("power")
    p.add_argument("--n", nargs="+", type=int, default=[5, 10, 15, 20, 25])
    p.add_argument("--cv", nargs="+", type=float, default=[0.2, 0.3])
    p.add_argument("--phi", nargs="+", type=float, default=[0.0, 0.3])
    p.add_argument("--sims", type=int, default=300); p.add_argument("--seed", type=int, default=0)
    p.set_defaults(fn=cmd_power)
    p = sub.add_parser("clusters")
    p.add_argument("--kmin", type=int, default=4); p.add_argument("--kmax", type=int, default=10)
    p.add_argument("--k", type=int, help="fijar k (si no, el de mayor silueta)")
    p.add_argument("--seed", type=int, default=0); p.set_defaults(fn=cmd_clusters)
    p = sub.add_parser("field"); p.set_defaults(fn=cmd_field)

    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()

"""Catálogo de medidas para análisis autoguiado ("una talla para todos").

En lugar de preguntas prefabricadas, el usuario elige:
    LUGAR      punto, área dibujada o entidad INEGI
    MEDIDA     una del CATALOG (abajo)
    MESES      cualquier subconjunto (temporada), p. ej. [6, 7, 8, 9]
    PERIODO    años a describir
    CONTRA     referencia climatológica, otro periodo u otro lugar

y cada resultado regresa con su incertidumbre, su prueba y un NIVEL DE CONFIANZA que depende
de cuántos años la sostienen (no de lo que el usuario quiera ver).

    python measures.py list
    python measures.py total --cve 09 --months 6 7 8 9 --period 2021 2025
    python measures.py drought --cve 09 --scale 3 --period 2020 2025
    python measures.py yoy --cve 14067 --months 7 --period 2024 2024
    python measures.py compare_periods --cve 09 --months 6 7 8 9 --a 1998 2010 --b 2011 2025

Fuentes: mensual GPM_3IMERGM (1998→) para medidas de escala mensual; diario GPM_3IMERGDF
(2020→) para medidas de escala diaria. No se mezclan en una misma serie (difieren ~3–4 %).
"""
import argparse
import json
import warnings
from dataclasses import asdict, dataclass, field
from functools import lru_cache

import numpy as np
import pandas as pd
import pymannkendall as mk
import xarray as xr
from scipy import stats

import config
from analysis import ALPHA, HEAVY_MM, MONTHS_ES, WET_MM, _json_default, spi_class

BASELINE = (1998, 2020)    # referencia por defecto (23 años de IMERG; OMM pide 30, 1991–2020)


# ---------------------------------------------------------------------------
# Catálogo
# ---------------------------------------------------------------------------
@dataclass
class Measure:
    id: str
    nombre: str
    pregunta: str
    unidad: str
    datos: str                 # "mensual" | "diario"
    prueba: str                # método de incertidumbre / significancia
    alta: int = 20             # años para confianza alta
    media: int = 8             # años para confianza media (debajo: baja, sólo descriptivo)
    parametros: list = field(default_factory=list)
    advertencia: str = ""


CATALOG = {m.id: m for m in [
    Measure("total", "Lluvia acumulada", "¿Cuánto llueve en estos meses?", "mm", "mensual",
            "Media y mediana entre años; IC95 bootstrap por años; error aleatorio IMERG como cota."),
    Measure("pct_normal", "% de lo normal", "¿Este periodo fue más seco o más húmedo de lo normal?", "%", "mensual",
            "Total / mediana de la referencia; percentil empírico y rango histórico (el más seco desde…).",
            advertencia="La referencia (1998–2020) tiene 23 años; percentiles extremos (<5 % o >95 %) son poco precisos."),
    Measure("drought", "Probabilidad de sequía (SPI)", "¿Qué tan seguido hay sequía meteorológica?", "% de meses", "mensual",
            "SPI-k (gamma por mes calendario con masa en cero, ajustada a la referencia); frecuencia de SPI ≤ −1 "
            "con IC de Wilson sobre n efectivo (corrige el traslape de ventanas) y binomial vs 15.9 % esperado.",
            parametros=["scale (1, 3, 6, 12 meses)"]),
    Measure("yoy", "Cambio interanual", "¿Cuánto cambió respecto al mismo periodo del año anterior?", "%", "mensual",
            "Cambio % año contra año; se compara con la distribución de cambios interanuales de la referencia "
            "(¿es un cambio inusual o variabilidad normal?)."),
    Measure("mom", "Cambio mes a mes vs lo normal", "¿El cambio respecto al mes anterior fue mayor de lo habitual?", "%",
            "mensual", "Cambio % observado vs el cambio climatológico normal entre esos dos meses; percentil del "
            "cociente en la referencia.",
            advertencia="El cambio mes a mes crudo está dominado por el ciclo estacional (junio siempre es +200 % "
                        "sobre mayo en el centro del país); por eso se reporta contra el cambio normal."),
    Measure("trend", "Tendencia de largo plazo", "¿Está aumentando o disminuyendo la lluvia?", "%/década", "mensual",
            "Pendiente de Sen + Hamed-Rao (autocorrelación) en totales anuales de los meses elegidos; MK estacional si "
            "son varios meses; prueba de Pettitt de cambio abrupto; tendencia mínima detectable si no es significativa.",
            alta=25, media=15,
            advertencia="IMERG cambia de constelación en 2014 (TRMM→GPM); si Pettitt marca ~2014, sospechar artefacto."),
    Measure("compare_periods", "Comparar dos periodos", "¿Cambió la lluvia entre el periodo A y el B?", "mm", "mensual",
            "Mann-Whitney (exacta) sobre totales anuales; desplazamiento de Hodges-Lehmann con IC95; delta de Cliff.",
            alta=10, media=5, parametros=["a (años)", "b (años)"]),
    Measure("compare_places", "Comparar dos lugares", "¿Llueve distinto en el lugar A que en el B?", "mm", "mensual",
            "Wilcoxon de rangos con signo sobre diferencias pareadas por año; mediana de la diferencia con IC95.",
            alta=15, media=8, parametros=["cve2 / punto2"]),
    Measure("wet_days", "Probabilidad de día con lluvia", "¿Qué tan probable es que llueva un día cualquiera?", "%",
            "diario", "Proporción de días ≥ umbral; IC95 bootstrap por bloques (o por años si hay ≥ 8).",
            alta=15, media=8, parametros=["umbral mm (1 = lluvia, 10 = fuerte, 20 = muy fuerte)"]),
    Measure("dry_spell", "Rachas secas", "¿Qué tan probable es una racha seca de N días?", "%", "diario",
            "Probabilidad de que una ventana de N días no tenga lluvia; racha seca más larga (CDD) por año.",
            alta=15, media=8, parametros=["dias (N)"]),
    Measure("intensity", "Intensidad y extremos", "¿Qué tan fuertes son los aguaceros?", "mm", "diario",
            "SDII (mm por día con lluvia), Rx1day y Rx5day por año: mediana y rango.", alta=15, media=8,
            advertencia="IMERG a 0.1° suaviza los extremos locales; Rx1day de una estación será mayor."),
]}


def confidence(n_years: int, m: Measure) -> str:
    if n_years < 3:
        return "no disponible"
    if n_years < m.media:
        return "baja"
    return "alta" if n_years >= m.alta else "media"


CONF_NOTE = {
    "alta": "Suficientes años para estimar la variabilidad entre años y probar significancia.",
    "media": "Resultado útil pero con IC amplios; la significancia sólo detecta diferencias grandes.",
    "baja": "Sólo descriptivo: muy pocos años para afirmar diferencias o tendencias (no se reporta significancia).",
    "no disponible": "No hay suficientes años con datos completos para esta medida.",
}


# ---------------------------------------------------------------------------
# Series por lugar
# ---------------------------------------------------------------------------
@lru_cache(maxsize=2)
def _open_cached(path: str, mtime: float) -> xr.Dataset:
    return xr.open_dataset(path).load()


def _monthly_grid() -> xr.Dataset:
    """Se recarga sola si download_monthly.py reescribe el archivo (clave = ruta + mtime)."""
    if not config.MONTHLY_FILE.exists():
        raise SystemExit("Falta imerg_mexico_monthly.nc: corre `python download_monthly.py`.")
    return _open_cached(str(config.MONTHLY_FILE), config.MONTHLY_FILE.stat().st_mtime)


def _monthly_zonal() -> xr.Dataset:
    return _open_cached(str(config.ZONAL_MONTHLY), config.ZONAL_MONTHLY.stat().st_mtime)


def _nearest_valid(da: xr.DataArray, lat: float, lon: float) -> tuple[int, int]:
    valid = da.notnull().any("time").values
    LA, LO = np.meshgrid(da.lat.values, da.lon.values, indexing="ij")
    d2 = np.where(valid, (LA - lat) ** 2 + ((LO - lon) * np.cos(np.deg2rad(lat))) ** 2, np.inf)
    return np.unravel_index(np.argmin(d2), d2.shape)


def _polygon_mean(da: xr.DataArray, geometry: dict) -> pd.Series:
    from shapely.geometry import shape
    from zonal import polygon_weights
    ii, jj, frac = polygon_weights(shape(geometry), da.lat.values, da.lon.values)
    if not len(ii):
        raise ValueError("El área está fuera de la malla de México.")
    w = frac * np.cos(np.deg2rad(da.lat.values[ii]))
    vals = da.values[:, ii, jj]
    ok = ~np.isnan(vals)
    with np.errstate(invalid="ignore"):
        s = (np.nan_to_num(vals) * w).sum(1) / (ok * w).sum(1)
    return pd.Series(s, index=da.time.to_index())


def monthly_series(place: dict, var: str = "precipitation") -> tuple[pd.Series, str]:
    """mm/mes para {'cve': '09'} | {'point': [lat, lon]} | {'geometry': GeoJSON}."""
    if "cve" in place:
        z = _monthly_zonal()
        r = z.sel(region=place["cve"])
        return r[var].to_series(), str(r.nombre.values)
    da = _monthly_grid()[var]
    if "point" in place:
        i, j = _nearest_valid(da, *place["point"])
        return da.isel(lat=i, lon=j).to_series(), f"Celda ({float(da.lat[i]):.2f}, {float(da.lon[j]):.2f})"
    return _polygon_mean(da, place["geometry"]), place.get("name", "Área")


_DAILY_SOURCE = None


def set_daily_source(fn) -> None:
    """Permite a la API inyectar su propio lector (con caché) de series diarias por lugar."""
    global _DAILY_SOURCE
    _DAILY_SOURCE = fn


def available(m: Measure) -> tuple[bool, str]:
    if m.datos == "mensual" and not (config.MONTHLY_FILE.exists() and config.ZONAL_MONTHLY.exists()):
        return False, "Faltan datos mensuales 1998→ (python download_monthly.py)"
    return True, ""


def daily_series(place: dict) -> tuple[pd.Series, str]:
    if _DAILY_SOURCE:
        return _DAILY_SOURCE(place)
    from analysis import open_grid, point_series, region_series
    if "cve" in place:
        return region_series([place["cve"]])
    if "point" in place:
        return point_series(*place["point"])
    pr = open_grid()["precipitation"].load()
    return _polygon_mean(pr, place["geometry"]), place.get("name", "Área")


def yearly_totals(m: pd.Series, months: list[int]) -> pd.Series:
    """Total de los `months` por año (sólo años con todos esos meses). Temporadas que cruzan el año
    (p. ej. nov–mar) se asignan al año en que terminan."""
    months = list(months)
    x = m[m.index.month.isin(months)].dropna()
    wraps = months != sorted(months)
    key = x.index.year + ((x.index.month >= months[0]) & wraps).astype(int) if wraps else x.index.year
    g = x.groupby(key)
    return g.sum()[g.size() == len(months)]


def _years(s: pd.Series, period) -> pd.Series:
    return s[(s.index >= period[0]) & (s.index <= period[1])] if period else s


def _boot_ci(x: np.ndarray, fn=np.mean, n=4000, seed=0) -> list:
    if len(x) < 3:
        return [None, None]
    rng = np.random.default_rng(seed)
    b = np.array([fn(rng.choice(x, len(x))) for _ in range(n)])
    return [round(float(np.percentile(b, 2.5)), 1), round(float(np.percentile(b, 97.5)), 1)]


def _p(p: float) -> float:
    """p-valor con 2 cifras significativas (0.0000012 en lugar de 0.0)."""
    return float(f"{p:.2g}")


def _result(mid: str, label: str, n_years: int, **kw) -> dict:
    m = CATALOG[mid]
    conf = confidence(n_years, m)
    out = {"medida": m.nombre, "lugar": label, "n_anios": n_years, "confianza": conf,
           "nota_confianza": CONF_NOTE[conf], "metodo": m.prueba, **kw}
    if m.advertencia:
        out["advertencia"] = m.advertencia
    if conf in ("baja", "no disponible"):   # no se reporta significancia sin potencia
        for k in [k for k in out if k.startswith("p_") or k == "significativo"]:
            out.pop(k)
    return out


# ---------------------------------------------------------------------------
# Medidas mensuales
# ---------------------------------------------------------------------------
def m_total(place, months=range(1, 13), period=None, **_):
    m, label = monthly_series(place)
    tot = _years(yearly_totals(m, months), period)
    err = None
    try:
        e, _ = monthly_series(place, "random_error")
        err = round(float(_years(yearly_totals(e, months), period).mean()), 1)
    except (KeyError, SystemExit):
        pass
    return _result("total", label, len(tot), meses=_mlabel(months),
                   anios=[int(tot.index.min()), int(tot.index.max())] if len(tot) else None,
                   media_mm=round(float(tot.mean()), 1), mediana_mm=round(float(tot.median()), 1),
                   ic95_media_mm=_boot_ci(tot.values), min_mm=round(float(tot.min()), 1),
                   max_mm=round(float(tot.max()), 1),
                   error_aleatorio_imerg_mm=err, por_anio=tot.round(1).to_dict())


def m_pct_normal(place, months=range(1, 13), period=None, baseline=BASELINE, **_):
    m, label = monthly_series(place)
    tot = yearly_totals(m, months)
    ref = _years(tot, baseline)
    ev = _years(tot, period)
    med = ref.median()
    rows = {}
    for y, v in ev.items():
        drier = tot[(tot.index < y) & (tot <= v)]
        wetter = tot[(tot.index < y) & (tot >= v)]
        rows[int(y)] = {"total_mm": round(float(v), 1), "pct_normal": round(float(100 * v / med), 1),
                        "percentil_ref": round(float(stats.percentileofscore(ref.values, v, kind="mean")), 1),
                        "clase": _tercile(v, ref.values),
                        # año más reciente igual o más extremo; "récord" si no hubo ninguno antes en el registro
                        "mas_seco_desde": (int(drier.index.max()) if len(drier) else f"récord desde {int(tot.index.min())}")
                        if v < med else None,
                        "mas_humedo_desde": (int(wetter.index.max()) if len(wetter) else f"récord desde {int(tot.index.min())}")
                        if v > med else None}
    return _result("pct_normal", label, len(ref), meses=_mlabel(months), referencia=list(baseline),
                   normal_mm=round(float(med), 1),
                   rango_normal_p33_p67=[round(float(np.percentile(ref, 33.3)), 1), round(float(np.percentile(ref, 66.7)), 1)],
                   por_anio=rows)


def _tercile(v, ref):
    lo, hi = np.percentile(ref, [33.3, 66.7])
    return "abajo de lo normal" if v < lo else "arriba de lo normal" if v > hi else "normal"


def spi_series(m: pd.Series, scale: int, baseline=BASELINE) -> pd.Series:
    """SPI-k: acumulado móvil de k meses; gamma por mes calendario ajustada a la referencia."""
    acc = m.rolling(scale, min_periods=scale).sum()
    out = pd.Series(np.nan, index=acc.index)
    for mo in range(1, 13):
        sel = acc.index.month == mo
        x = acc[sel]
        ref = x[(x.index.year >= baseline[0]) & (x.index.year <= baseline[1])].dropna().values
        if len(ref) < 10:
            continue
        q = np.mean(ref <= 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a, _, sc = stats.gamma.fit(ref[ref > 0], floc=0)
        cdf = q + (1 - q) * stats.gamma.cdf(np.clip(x.values, 0, None), a, scale=sc)
        out[sel] = stats.norm.ppf(np.clip(cdf, 1e-6, 1 - 1e-6))
    return out


def m_drought(place, period=None, baseline=BASELINE, scale=3, **_):
    m, label = monthly_series(place)
    s = spi_series(m, scale, baseline).dropna()
    ref_years = s[(s.index.year >= baseline[0]) & (s.index.year <= baseline[1])].index.year.nunique()
    ev = s[(s.index.year >= period[0]) & (s.index.year <= period[1])] if period else s
    k, n = int((ev <= -1).sum()), len(ev)
    r1 = max(0.0, ev.autocorr(1)) if n > 3 else 0.0
    n_eff = max(1.0, n * (1 - r1) / (1 + r1))                       # meses efectivamente independientes
    p_hat = k / n if n else np.nan
    z = 1.96   # IC de Wilson con n efectivo
    center = (p_hat + z * z / (2 * n_eff)) / (1 + z * z / n_eff) if n else np.nan
    half = z * np.sqrt(p_hat * (1 - p_hat) / n_eff + z * z / (4 * n_eff ** 2)) / (1 + z * z / n_eff) if n else np.nan
    k_eff = p_hat * n_eff
    runs = (ev <= -1).astype(int)
    longest = int(runs.groupby((runs == 0).cumsum()).sum().max()) if n else 0
    p0 = stats.norm.cdf(-1)
    p_val = stats.binomtest(int(round(k_eff)), int(round(n_eff)), p0, "greater").pvalue if n else None
    last = ev.iloc[-1] if n else np.nan
    return _result("drought", label, int(ref_years), escala_meses=scale, referencia=list(baseline),
                   meses_evaluados=n, n_efectivo=round(n_eff, 1),
                   prob_sequia_pct=round(100 * p_hat, 1),
                   ic95_pct=[round(100 * max(0, center - half), 1), round(100 * min(1, center + half), 1)],
                   esperado_pct=round(100 * p0, 1),
                   p_mas_que_esperado=None if p_val is None else round(float(p_val), 4),
                   racha_sequia_mas_larga_meses=longest,
                   ultimo_mes={"fecha": str(ev.index[-1].date()) if n else None,
                               "SPI": None if np.isnan(last) else round(float(last), 2),
                               "clase": None if np.isnan(last) else str(spi_class(np.array([last]))[0])})


def m_yoy(place, months=range(1, 13), period=None, baseline=BASELINE, **_):
    m, label = monthly_series(place)
    tot = yearly_totals(m, months)
    ref = _years(tot, baseline)
    ref_chg = (ref / ref.shift(1) - 1).dropna() * 100                   # cambios interanuales "normales"
    rows = {}
    for y in _years(tot, period).index:
        if y - 1 not in tot.index:
            continue
        chg = 100 * (tot[y] / tot[y - 1] - 1)
        rows[int(y)] = {"total_mm": round(float(tot[y]), 1), "anterior_mm": round(float(tot[y - 1]), 1),
                        "cambio_pct": round(float(chg), 1),
                        "percentil_vs_cambios_normales": round(float(stats.percentileofscore(ref_chg, chg, kind="mean")), 1),
                        "inusual": bool(chg < np.percentile(ref_chg, 10) or chg > np.percentile(ref_chg, 90))}
    return _result("yoy", label, len(ref_chg), meses=_mlabel(months), referencia=list(baseline),
                   cambio_normal_p10_p90_pct=[round(float(np.percentile(ref_chg, 10)), 1),
                                              round(float(np.percentile(ref_chg, 90)), 1)] if len(ref_chg) else None,
                   por_anio=rows)


def m_mom(place, period=None, baseline=BASELINE, **_):
    m, label = monthly_series(place)
    m = m.dropna()
    ratio = (m / m.shift(1)).where(m.shift(1) > 0)
    ref = ratio[(ratio.index.year >= baseline[0]) & (ratio.index.year <= baseline[1])]
    normal = ref.groupby(ref.index.month).median()                     # cambio típico entre esos dos meses
    ev = ratio[(ratio.index.year >= period[0]) & (ratio.index.year <= period[1])] if period else ratio.iloc[-12:]
    rows = {}
    for t, r in ev.dropna().items():
        refm = ref[ref.index.month == t.month].dropna()
        rows[str(t.date())[:7]] = {
            "mm": round(float(m[t]), 1), "cambio_pct": round(float(100 * (r - 1)), 1),
            "cambio_normal_pct": round(float(100 * (normal[t.month] - 1)), 1),
            "exceso_vs_normal_pct": round(float(100 * (r / normal[t.month] - 1)), 1),
            "percentil": round(float(stats.percentileofscore(np.log(refm), np.log(r), kind="mean")), 1)}
    return _result("mom", label, int(ref.index.year.nunique()), referencia=list(baseline), por_mes=rows)


def pettitt(x: np.ndarray) -> tuple[int, float]:
    n = len(x)
    sgn = np.sign(x[None, :] - x[:, None])
    U = np.array([sgn[:t + 1, t + 1:].sum() for t in range(n - 1)])
    K = np.abs(U).max()
    return int(np.abs(U).argmax()), float(min(1.0, 2 * np.exp(-6 * K ** 2 / (n ** 3 + n ** 2))))


def m_trend(place, months=range(1, 13), period=None, **_):
    from analysis import mde
    m, label = monthly_series(place)
    tot = _years(yearly_totals(m, months), period)
    n = len(tot)
    out = {"meses": _mlabel(months), "anios": [int(tot.index.min()), int(tot.index.max())] if n else None}
    if n >= 8:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hr = mk.hamed_rao_modification_test(tot.values, alpha=ALPHA)
        mean = tot.mean()
        out.update(sen_mm_por_decada=round(10 * hr.slope, 1), sen_pct_por_decada=round(1000 * hr.slope / mean, 1),
                   p_hamed_rao=_p(hr.p), significativo=bool(hr.p < ALPHA), direccion=hr.trend)
        if len(list(months)) > 1:
            mm = m[m.index.month.isin(list(months))]
            mm = mm[mm.index.year.isin(tot.index)]
            st = mk.seasonal_test(mm.values, period=len(list(months)), alpha=ALPHA)
            out["p_mk_estacional"] = round(float(st.p), 4)
        idx, pp = pettitt(tot.values)
        out["pettitt"] = {"anio_cambio": int(tot.index[idx + 1]), "p": round(pp, 4),
                          "nota": "cambio abrupto cerca de 2014 puede ser artefacto TRMM→GPM" if pp < ALPHA else None}
        if hr.p >= ALPHA:
            cv = float(tot.std() / mean)
            b = mde(n, round(cv, 2), max(0.0, float(tot.autocorr(1))), 200, np.random.default_rng(0))
            out["tendencia_minima_detectable_pct_decada"] = None if np.isinf(b) else round(1000 * b, 0)
            out["lectura"] = ("Sin tendencia significativa. Esto NO prueba que no haya cambio: tendencias menores a "
                              f"~{out['tendencia_minima_detectable_pct_decada']} %/década no se pueden detectar con {n} años.")
    return _result("trend", label, n, **out)


def m_compare_periods(place, months=range(1, 13), a=None, b=None, **_):
    m, label = monthly_series(place)
    tot = yearly_totals(m, months)
    A, B = _years(tot, a).values, _years(tot, b).values
    n = min(len(A), len(B))
    out = {"meses": _mlabel(months), "a": list(a), "b": list(b), "n_a": len(A), "n_b": len(B),
           "mediana_a_mm": round(float(np.median(A)), 1), "mediana_b_mm": round(float(np.median(B)), 1)}
    if n >= 3:
        diffs = (B[:, None] - A[None, :]).ravel()
        hl = float(np.median(diffs))                                    # Hodges-Lehmann
        mw = stats.mannwhitneyu(B, A, alternative="two-sided", method="auto")
        cliff = float((np.sign(diffs)).mean())
        rng = np.random.default_rng(0)
        boot = [np.median((rng.choice(B, len(B))[:, None] - rng.choice(A, len(A))[None, :]).ravel()) for _ in range(2000)]
        out.update(cambio_b_menos_a_mm=round(hl, 1), cambio_pct=round(float(100 * hl / np.median(A)), 1),
                   ic95_mm=[round(float(np.percentile(boot, 2.5)), 1), round(float(np.percentile(boot, 97.5)), 1)],
                   p_mann_whitney=_p(mw.pvalue), significativo=bool(mw.pvalue < ALPHA),
                   delta_cliff=round(cliff, 2),
                   tamano_efecto="despreciable" if abs(cliff) < 0.147 else "pequeño" if abs(cliff) < 0.33
                   else "mediano" if abs(cliff) < 0.474 else "grande")
    return _result("compare_periods", label, n, **out)


def m_compare_places(place, place2=None, months=range(1, 13), period=None, **_):
    m1, l1 = monthly_series(place)
    m2, l2 = monthly_series(place2)
    t1, t2 = _years(yearly_totals(m1, months), period), _years(yearly_totals(m2, months), period)
    common = t1.index.intersection(t2.index)
    d = (t2[common] - t1[common]).values
    out = {"lugar_b": l2, "meses": _mlabel(months), "mediana_a_mm": round(float(t1[common].median()), 1),
           "mediana_b_mm": round(float(t2[common].median()), 1)}
    if len(d) >= 5:
        w = stats.wilcoxon(d)
        out.update(diferencia_mediana_b_menos_a_mm=round(float(np.median(d)), 1), ic95_mm=_boot_ci(d, np.median),
                   p_wilcoxon=_p(w.pvalue), significativo=bool(w.pvalue < ALPHA),
                   anios_b_mas_lluvioso=int((d > 0).sum()), anios_comparados=len(d))
    return _result("compare_places", l1, len(d), **out)


# ---------------------------------------------------------------------------
# Medidas diarias
# ---------------------------------------------------------------------------
def _daily_by_year(place, months, period):
    s, label = daily_series(place)
    s = s[s.index.month.isin(list(months))]
    if period:
        s = s[(s.index.year >= period[0]) & (s.index.year <= period[1])]
    days_per_year = s.notna().groupby(s.index.year).sum()
    full = s.index.to_series().groupby(s.index.year).count()
    good = days_per_year[days_per_year >= 0.9 * full].index
    return s[s.index.year.isin(good)], label, list(good)


def m_wet_days(place, months=range(1, 13), period=None, threshold=WET_MM, **_):
    from analysis import MIN_YEARS_YEAR_BOOT, bootstrap_stat
    s, label, years = _daily_by_year(place, months, period)
    ind = (s >= threshold).astype(float).where(s.notna())
    segs = [g.dropna().values for _, g in ind.groupby(ind.index.year)]
    by_year = len(years) >= MIN_YEARS_YEAR_BOOT
    b = bootstrap_stat(segs, np.mean, 7, 2000, np.random.default_rng(0), by_year) if segs else np.array([np.nan])
    p = float(ind.mean())
    return _result("wet_days", label, len(years), meses=_mlabel(months), umbral_mm=threshold,
                   prob_pct=round(100 * p, 1), ic95_pct=[round(100 * np.nanpercentile(b, 2.5), 1),
                                                         round(100 * np.nanpercentile(b, 97.5), 1)],
                   dias_esperados_por_mes=round(30.4 * p, 1), anios=years,
                   ic_incluye_variabilidad_interanual=by_year)


def m_dry_spell(place, months=range(1, 13), period=None, dias=7, **_):
    s, label, years = _daily_by_year(place, months, period)
    wet = (s >= WET_MM).astype(float).where(s.notna())
    win = wet[::-1].rolling(dias, min_periods=dias).sum()[::-1]
    p = float((win == 0).mean())
    cdd = {}
    for y, g in wet.groupby(wet.index.year):
        dry = (g == 0).astype(int)
        cdd[int(y)] = int(dry.groupby((dry == 0).cumsum()).sum().max())
    return _result("dry_spell", label, len(years), meses=_mlabel(months), dias=dias,
                   prob_ventana_seca_pct=round(100 * p, 1),
                   racha_seca_mas_larga_por_anio=cdd, racha_mediana=float(np.median(list(cdd.values()))) if cdd else None)


def m_intensity(place, months=range(1, 13), period=None, **_):
    s, label, years = _daily_by_year(place, months, period)
    rows = {}
    for y, g in s.groupby(s.index.year):
        rows[int(y)] = {"SDII_mm": round(float(g[g >= WET_MM].mean()), 1), "Rx1day_mm": round(float(g.max()), 1),
                        "Rx5day_mm": round(float(g.rolling(5).sum().max()), 1),
                        "dias_fuertes": int((g >= HEAVY_MM).sum())}
    df = pd.DataFrame(rows).T
    return _result("intensity", label, len(years), meses=_mlabel(months),
                   mediana=df.median().round(1).to_dict() if len(df) else None, por_anio=rows)


FUNCS = {"total": m_total, "pct_normal": m_pct_normal, "drought": m_drought, "yoy": m_yoy, "mom": m_mom,
         "trend": m_trend, "compare_periods": m_compare_periods, "compare_places": m_compare_places,
         "wet_days": m_wet_days, "dry_spell": m_dry_spell, "intensity": m_intensity}


CELL_DEG = 0.1


def resolution_note(place: dict) -> dict | None:
    """Advierte cuando el polígono es del orden de una celda IMERG (0.1° ≈ 11 km): el resultado
    representa a la(s) celda(s), no al detalle del polígono."""
    if "geometry" not in place:
        return None
    from shapely.geometry import shape
    from zonal import KM_PER_DEG, polygon_weights
    g = shape(place["geometry"])
    lat0 = g.centroid.y
    area_km2 = g.area * KM_PER_DEG ** 2 * np.cos(np.deg2rad(lat0))
    cell_km2 = (CELL_DEG * KM_PER_DEG) ** 2 * np.cos(np.deg2rad(lat0))
    grid = _monthly_grid() if config.MONTHLY_FILE.exists() else None
    lat = grid.lat.values if grid is not None else np.arange(14.45, 32.85, CELL_DEG)
    lon = grid.lon.values if grid is not None else np.arange(-118.45, -86.55, CELL_DEG)
    ii, _, frac = polygon_weights(g, lat, lon)
    note = {"area_km2": round(float(area_km2), 1), "celdas_tocadas": int(len(ii)),
            "area_equivalente_celdas": round(float(area_km2 / cell_km2), 2)}
    if area_km2 < 3 * cell_km2:
        size = ("menos de una celda" if area_km2 < cell_km2 else f"unas {area_km2 / cell_km2:.1f} celdas")
        note["aviso"] = (f"El área ({area_km2:,.0f} km²) equivale a {size} IMERG de ~{cell_km2:,.0f} km² "
                         f"(toca {len(ii)}). El resultado describe esas celdas completas, no variaciones dentro del polígono.")
    return note


def compute(measure: str, place: dict, **params) -> dict:
    res = FUNCS[measure](place, **params)
    note = resolution_note(place)
    if note:
        res["resolucion"] = note
    return res


def _mlabel(months) -> str:
    months = list(months)
    return "año completo" if len(months) == 12 else "–".join(MONTHS_ES[m - 1] for m in (months[0], months[-1])) \
        if len(months) > 2 else ", ".join(MONTHS_ES[m - 1] for m in months)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("measure", choices=["list", *FUNCS])
    ap.add_argument("--cve"); ap.add_argument("--point", nargs=2, type=float)
    ap.add_argument("--cve2"); ap.add_argument("--point2", nargs=2, type=float)
    ap.add_argument("--months", nargs="+", type=int, default=list(range(1, 13)))
    ap.add_argument("--period", nargs=2, type=int); ap.add_argument("--baseline", nargs=2, type=int, default=BASELINE)
    ap.add_argument("--a", nargs=2, type=int); ap.add_argument("--b", nargs=2, type=int)
    ap.add_argument("--scale", type=int, default=3); ap.add_argument("--threshold", type=float, default=WET_MM)
    ap.add_argument("--dias", type=int, default=7)
    a = ap.parse_args()
    if a.measure == "list":
        for m in CATALOG.values():
            print(f"{m.id:16s} [{m.datos}] {m.nombre} — {m.pregunta}  (alta ≥{m.alta} años, media ≥{m.media})")
        return
    place = {"cve": a.cve} if a.cve else {"point": a.point}
    place2 = {"cve": a.cve2} if a.cve2 else ({"point": a.point2} if a.point2 else None)
    res = compute(a.measure, place, months=a.months, period=a.period, baseline=tuple(a.baseline), a=a.a, b=a.b,
                  scale=a.scale, threshold=a.threshold, dias=a.dias, place2=place2)
    print(json.dumps(res, ensure_ascii=False, indent=2, default=_json_default))


if __name__ == "__main__":
    main()

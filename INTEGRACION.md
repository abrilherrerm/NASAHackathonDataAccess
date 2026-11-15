# Integración al repositorio final

Guía para llevar el trabajo de este repositorio (interino) al repo final **sin tocar la interfaz nueva**
ni sobrescribir cambios del equipo. Todo lo descrito está en `dev` (commit `e6103de` o posterior).

**Regla de oro:** no copiar `app/static/index.html` de aquí. Todo lo demás es archivo nuevo o un parche
pequeño y localizado.

---

## 1. Copiar tal cual (archivos nuevos, no chocan con nada)

| Archivo | Qué es |
|---|---|
| `download_monthly.py` | Descarga y proceso del mensual IMERG (GPM_3IMERGM) 1998 → hoy |
| `data/processed/imerg_mexico_monthly.nc` (41 MB) | Mensual 1998-01 → 2025-09: precipitación, error aleatorio, calidad |
| `data/processed/zonal_monthly.nc` (7 MB) | Series mensuales de México, 32 estados y 2,478 municipios |
| `analysis.py` | Pruebas estadísticas (forma de los datos, Mann-Kendall/Hamed-Rao/Sen, SPI, FDR, regiones) |
| `measures.py` | Catálogo de 11 medidas con nivel de confianza (lo que usa la API) |
| `validate_late.py` | Validación IMERG Late vs Final (justifica excluir Late de la estadística) |
| `ANALISIS.md` | Métodos, decisiones y evidencia (§8 temporadas, §9 Late) |
| `data/analysis/*` | Resultados regenerados con los datos completos |

## 2. Añadir a archivos existentes (líneas al final, sin borrar nada)

- **`config.py`** — bloque al final (`MONTHLY_CONCEPT_ID`, `MONTHLY_START`, `RAW_MONTHLY_DIR`, `MONTHLY_FILE`, `ZONAL_MONTHLY`).
- **`requirements.txt`** — `pymannkendall>=1.4.3`.
- **`.gitignore`** — `data/raw_monthly` y `data/raw_late_overlap`.

## 3. Parche a `app/main.py` (revisar contra la versión del repo final)

**Agregar (no cambia nada existente):**
- Endpoints `GET /api/measures` y `POST /api/measure`, con `MeasureRequest`, `_daily_for_place`, `_clean`
  y la línea `measures.set_daily_source(_daily_for_place)`.
- `_valid_cells` y `_nearest_cell`: un clic en la costa usa la celda con dato más cercana (antes daba 404).

**Cambia comportamiento — decidir si se adopta:**
- `_slice_final` en los 4 endpoints de series (`/api/timeseries/point`, `region`, `polygon`, `entity/{cvegeo}`):
  las series terminan en el último día **Final** (hoy 2025-09-30). `/api/map` **no** cambia: los mapas siguen
  mostrando IMERG Late (preliminar). Motivo: Late sobreestima ~+18 % regional y no se puede corregir (ANALISIS.md §9).
- `/api/info` agrega `final_end`; `_late_from` devuelve `null` cuando las series ya están recortadas.
- `point_series` usa `_nearest_cell` en lugar de `sel(..., method="nearest")`.

## 4. Interfaz: funciones a reconectar en la UI nueva (no copiar HTML)

La UI interina tenía estas funciones; la UI nueva puede implementarlas como quiera usando la API:

1. **Comparación A | B** — dos selecciones lado a lado; con escala compartida (misma escala de color en mapas,
   mismo eje Y en series) y diferencias de B vs A.
2. **Tarjeta "Medidas"** — medida + temporada/meses + años, sobre la selección actual (punto, área, entidad).
3. **Temporadas** — presets "lluvias" / "secas" / año completo / núcleo jun–sep / personalizado.
4. **Fechas de series** con máximo = `info.final_end`; mapas con máximo = `info.end`.

### Contrato de la API de medidas

`GET /api/measures` → lista del catálogo:
`id, nombre, pregunta, unidad, datos ("diario"|"mensual"), prueba, alta, media, parametros, advertencia,
disponible, motivo, anios_disponibles [inicio, fin]`.

`POST /api/measure`:

```json
{
  "measure": "wet_days",
  "place": {"cve": "09"} | {"point": [19.43, -99.13]} | {"geometry": <GeoJSON>, "name": "Mi área"},
  "season": "lluvias" | "secas" | null,
  "months": [6, 7, 8, 9],
  "period": [2014, 2025],
  "threshold": 1, "dias": 7, "scale": 3,
  "a": [1998, 2011], "b": [2012, 2025], "place2": {"cve": "14067"}
}
```

- `season` manda sobre `months`. Diario: corte oficial exacto (15 may–15 nov / 16 nov–14 may).
  Mensual: meses completos (may–nov / dic–abr). La respuesta trae el nombre real en `meses`.
- Medidas diarias: `wet_days` (umbral), `dry_spell` (días), `intensity`.
  Mensuales: `total`, `pct_normal`, `drought` (escala SPI), `yoy`, `mom`, `trend`,
  `compare_periods` (a, b), `compare_places` (place2).

Campos comunes de la respuesta:

| Campo | Uso en la UI |
|---|---|
| `confianza` | `alta` / `media` / `baja` / `no disponible` → insignia |
| `nota_confianza` | texto que explica la insignia |
| `n_anios` | años completos que sostienen el resultado |
| `p_*`, `significativo` | **sólo vienen con confianza media/alta** (con baja se omiten a propósito) |
| `ic95_*` | intervalos de confianza |
| `por_anio` / `por_mes` | tabla de detalle |
| `advertencia`, `lectura` | avisos a mostrar tal cual |
| `resolucion.aviso` | polígonos menores a ~3 celdas IMERG: el valor describe celdas completas |
| `datos_hasta` | medidas diarias: último día Final usado |

**Convención para comparar A vs B:** cantidades (mm, días) en **% relativo**; probabilidades y porcentajes en
**puntos porcentuales (pp)**. No mostrar diferencia en el error IMERG ni en p-valores.

## 5. Verificar después de integrar

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload
curl -s localhost:8000/api/measures | head -c 300
curl -s -X POST localhost:8000/api/measure -H 'Content-Type: application/json' \
  -d '{"measure":"total","place":{"cve":"09"},"season":"lluvias"}'
# Esperado: "confianza": "alta", "n_anios": 27, "meses": "Temporada de lluvias (may–nov)"
python analysis.py travel --region pv     # Esperado: abril mejor mes; marzo y mayo indistinguibles
```

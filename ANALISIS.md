# Análisis estadístico — métodos y candados

`analysis.py` responde preguntas de precipitación con pruebas adecuadas a cómo es la lluvia
(muchos ceros, sesgo fuerte, estacionalidad y persistencia), y **se niega a reportar** lo que
los datos todavía no pueden sostener. Resultados en `data/analysis/`.

```bash
python analysis.py shape    --region pv        # forma de los datos -> qué métodos aplican
python analysis.py travel   --region pv        # ¿cuándo viajar a Puerto Vallarta?
python analysis.py season   --region zmvm      # temporadas de lluvia vs umbral de anomalía
python analysis.py trends   --region zmvm      # MK estacional / Hamed-Rao / Sen
python analysis.py power                       # tendencia mínima detectable por nº de años
python analysis.py clusters                    # regiones homogéneas + alto gradiente
python analysis.py field                       # tendencia por celda con control FDR
```
Objetivo: `--cve 14067` (claves INEGI), `--region pv|pv_bahia|zmvm` o `--point LAT LON`.

## 1. La forma de los datos decide (`shape`)

| Rasgo observado (PV, CDMX) | Consecuencia |
|---|---|
| 42–68 % de días secos, asimetría 2.8–4.2 | Nada de t-test ni OLS sobre mm. Pruebas de rangos (MK, Kruskal-Wallis), bootstrap, o gamma con masa en cero. |
| Totales mensuales: gamma gana por AIC | Respalda el **SPI** (gamma) para anomalías. Montos diarios: lognormal. |
| Estacionalidad (Kruskal-Wallis p < 0.01) | Comparar siempre contra el **mismo mes**; MK **estacional** para series mensuales. |
| Persistencia diaria r1 = 0.33–0.57 | Días no independientes: **bootstrap por bloques** (sino los IC salen demasiado estrechos → tipo I). |
| Autocorrelación interanual (si aparece) | **Hamed-Rao** en lugar de MK original (que infla el tipo I con r1 > 0). |

## 2. Candados contra errores tipo I y II

| Candado | Valor | Por qué |
|---|---|---|
| Día con lluvia | ≥ 1 mm | Umbral OMM/ETCCDI; IMERG produce llovizna espuria < 1 mm. |
| α por prueba | 0.05 | — |
| Comparaciones múltiples (12 meses) | **Holm** | Controla la tasa de error familiar sin ser tan conservador como Bonferroni. |
| Campos (miles de celdas) | **Benjamini-Hochberg, q = 0.10** | Con 17 mil celdas, 5 % salen "significativas" por azar. FDR robusta a correlación espacial (Wilks 2016). |
| Tendencias | **≥ 10 años** completos | Por debajo no hay potencia (ver tabla). |
| Climatología de referencia | **≥ 20 años** (OMM: 30) | Con menos, el umbral de anomalía se calcula con los mismos años que se evalúan. `--force` lo calcula marcado PROVISIONAL. |
| Bootstrap por años | **≥ 8 años** por mes | Sólo así los IC incluyen variabilidad interanual (ENSO, huracanes); antes se avisa que son optimistas. |
| Relevancia práctica | **h de Cohen** | Significativo ≠ importante: 0.2 pequeño · 0.5 mediano · 0.8 grande. |

**Potencia (`power`)**: tendencia mínima detectable con 80 % de potencia (Hamed-Rao, α = 0.05) en totales anuales gamma:

| CV anual | AR(1) | 5 años | 10 | 15 | 20 | 25 |
|---|---|---|---|---|---|---|
| 0.2 | 0.0 | no detectable | 78 %/déc | 41 | 25 | 18 |
| 0.3 | 0.3 | no detectable | 111 %/déc | 59 | 39 | 27 |

Con los 5 años actuales **no se puede afirmar ninguna tendencia**; ni con 10 se detectaría nada
menor a "la lluvia casi se duplicó". Con 25 años (IMERG desde 2000) se detectan cambios de ~20–27 %/década.

## 3. Preguntas y métricas mínimas viables

### ¿Cuándo viajar a Puerto Vallarta? (`travel`)
- **Métrica principal**: probabilidad de un *viaje seco* = ventana de N días (def. 7) con ≤ 1 día de lluvia.
  Secundarias: prob. diaria de lluvia (≥ 1 mm), de lluvia fuerte (≥ 10 mm), días de lluvia esperados.
- **Significancia**: bootstrap de la diferencia contra el mejor mes + Holm. Los meses *indistinguibles
  del mejor* forman la ventana recomendada; h de Cohen para ver si la diferencia importa.
- Resultado preliminar (celda IMERG, 2020-10 → 2025-09): **marzo–abril** ~100 % de viajes secos;
  febrero y mayo indistinguibles; nov–ene ~82–84 % (estadísticamente peor, h ≈ 0.8). Jul–sep ≈ 0 %.
- Pendiente: sólo precipitación. Temperatura/humedad y temporada de huracanes (jun–nov) cambiarían el ranking.

### ¿Cuántas temporadas de lluvia en la ZMVM rebasaron el umbral de anomalía? (`season`)
- **Región**: CDMX + 17 municipios conurbados del Edomex, promedio ponderado por área
  (aproximación; la delimitación oficial tiene 76 municipios — editar `REGIONS["zmvm"]`).
- **Temporada**: oficial 15 may–15 nov por defecto (`--window lluvias|secas`; `--months` para otra). Ver §8.
- **Umbral**: SPI de la temporada con gamma ajustada a la referencia; rebase = SPI ≤ −1 (seco) o ≥ +1 (húmedo),
  esperado 15.9 % de los años cada uno. También Rx1day > P90 de la referencia.
- **Frecuencia**: k/n con IC exacto de Clopper-Pearson + prueba binomial contra lo esperado.
- Índices ETCCDI por temporada: total, días de lluvia, SDII, Rx1day, Rx5day, R20mm, R95p, CDD.
- **Hoy no se puede responder**: no hay referencia 2001–2020 descargada.

### Tendencias (`trends`, `field`)
- Mensual: **MK estacional** (Hirsch 1982) y su versión correlacionada (Hirsch & Slack 1984).
- Anual y de temporada: **Hamed-Rao** (corrige autocorrelación) + MK original como contraste;
  si discrepan, manda Hamed-Rao. Pendiente de **Sen** (robusta a extremos).
- Campo: Hamed-Rao por celda + **FDR** → `field_trends.nc`.

## 4. Regiones homogéneas y alto gradiente (`clusters`)
- K-means sobre la **climatología** de cada celda: forma del ciclo anual (12 fracciones mensuales),
  log del total anual, frecuencia de días con lluvia, SDII e índice de estacionalidad. k por silueta (hoy k = 5).
- Se agrupa por **régimen**, no por tendencia: con pocos años, agrupar tendencias sería agrupar ruido.
  Después, la tendencia se prueba sobre la serie media de cada grupo (más potencia, menos pruebas).
- **Alto gradiente** = decil superior de |∇ log(total anual)| (montaña-costa, p. ej. Sierra Madre / Bahía de Banderas),
  y **frontera** = vecinos en otro grupo. Ahí la media del grupo no representa al sitio: usar celda o municipio.
- Silueta ~0.29: estructura moderada (los regímenes se mezclan gradualmente); tratar los bordes con cautela.

## 5. Qué falta para que las pruebas tengan potencia
1. ~~Ampliar a ~25 años~~ **Hecho para escala mensual**: `download_monthly.py` trae 1998-01 → 2025-09 (§6),
   suficiente para referencia, SPI y tendencias. Las medidas **diarias** siguen con pocos años: `config.START_DATE`
   ya apunta a 2014-06 (era GPM); al descargarlo pasan a confianza media (≥ 8 años) / alta (≥ 15).
2. Regenerar `zonal_series.nc` (`python zonal.py`) tras cada descarga del diario.
3. Validar IMERG contra estaciones del SMN en sitios de alto gradiente (IMERG subestima la lluvia orográfica).

---

## 6. Datos mensuales 1998 → 2025 (`download_monthly.py`)

**GPM_3IMERGM v07** (`C2723754851-GES_DISC`), misma malla de 0.1° que el diario. ~5 min de Harmony (1 job por año).

| Desbloquea (escala mensual) | No desbloquea (escala diaria) |
|---|---|
| Referencia 1998–2020 (23 años) para % de lo normal, terciles y SPI | Probabilidad de día con lluvia, rachas secas, Rx1day, "viaje seco" |
| SPI-1/3/6/12 (el SPI se define sobre totales mensuales) | Para estas: extender el **diario** a 2001 (~6 h de Harmony) |
| Tendencias con 27 años (detectan ~20–27 %/década) | |
| Climatología robusta para `clusters` | |

- **Unidades**: el archivo trae `precipitation` en **mm/h** (tasa media); se convierte a mm/mes × 24 × días del mes.
- **Consistencia**: vs la suma de los diarios (2024): r = 0.998 por celda-mes, el mensual ~3.6 % más alto
  (se calibra directo con pluviómetros GPCC). **No mezclar** productos en una serie: medidas mensuales usan el mensual.
- **Extras**: `random_error` (error aleatorio de IMERG, mm/mes) y `quality_index` / `gauge_weight`
  (cuántos pluviómetros respaldan la celda) — sirven para calificar la confianza espacial.
- **Homogeneidad**: IMERG pasa de TRMM a GPM en 2014; `trend` corre Pettitt para detectar saltos.

## 7. Medidas para análisis autoguiado (`measures.py`)

El usuario no elige una pregunta; elige **lugar × medida × meses × periodo × contra qué**. Cada resultado
trae valor, incertidumbre, prueba y un **nivel de confianza** fijado por los años que lo sostienen:

| Confianza | Regla | Qué se muestra |
|---|---|---|
| alta | ≥ `alta` años (20 en general; 25 tendencia) | valor + IC + significancia |
| media | ≥ `media` años (8; 15 tendencia) | valor + IC amplio + significancia (sólo detecta cambios grandes) |
| baja | < `media` | **sólo descriptivo**: se ocultan p-valores |
| no disponible | < 3 años | nada |

| Medida | Responde | Escala | Método | Hoy |
|---|---|---|---|---|
| `total` | ¿Cuánto llueve en estos meses? | mensual | media/mediana entre años, IC bootstrap por años, error IMERG | alta |
| `pct_normal` | ¿Más seco o húmedo de lo normal? | mensual | % de la mediana 1998–2020, percentil, terciles, "el más seco desde…" | alta |
| `drought` | ¿Qué tan seguido hay sequía? | mensual | SPI-k ≤ −1, IC de Wilson con n efectivo, binomial vs 15.9 % | alta |
| `yoy` | ¿Cuánto cambió vs el año pasado? | mensual | cambio %, comparado con los cambios interanuales normales | alta |
| `mom` | ¿El cambio mes a mes fue inusual? | mensual | cambio % **contra el cambio normal** entre esos meses | alta |
| `trend` | ¿Aumenta o disminuye? | mensual | Sen + Hamed-Rao, MK estacional, Pettitt, tendencia mínima detectable | alta |
| `compare_periods` | ¿Cambió entre A y B? | mensual | Mann-Whitney, Hodges-Lehmann + IC, delta de Cliff | alta |
| `compare_places` | ¿Llueve distinto en A y B? | mensual | Wilcoxon pareado por año, IC de la mediana | alta |
| `wet_days` | ¿Probabilidad de que llueva un día? | diario | proporción ≥ umbral (1/10/20 mm), bootstrap por bloques/años | baja |
| `dry_spell` | ¿Probabilidad de N días secos seguidos? | diario | ventanas sin lluvia, CDD por año | baja |
| `intensity` | ¿Qué tan fuertes son los aguaceros? | diario | SDII, Rx1day, Rx5day | baja |

Decisiones de diseño:
- **"% de aumento mes a mes" crudo no se ofrece**: el ciclo estacional lo domina (junio ≈ +200 % sobre mayo cada año).
  Se reporta el cambio observado **contra el cambio normal** y su percentil.
- **No significativo ≠ sin cambio**: `trend` reporta la tendencia mínima detectable para que la ausencia de
  evidencia no se lea como evidencia de ausencia (error tipo II).
- **Ocultar p-valores con confianza baja** evita que alguien encuentre "significancia" en 3 años de datos.
- Lugares: `{"cve": "09"}`, `{"point": [lat, lon]}` o `{"geometry": GeoJSON}` — los mismos que ya usa la app.

## 8. Temporadas: binario oficial lluvias / secas

**Por defecto se usa la temporada oficial**: lluvias **15 may–15 nov** y secas **16 nov–14 may** (complemento
exacto, sin traslape). Jun–sep queda como opción "núcleo lluvioso" sólo para comparar con la literatura del monzón.

**Por qué no jun–sep** (diario 2021–2024, `data/analysis/ventana_temporada.csv`):

| | Jun–sep | Oficial | Lluvia de temporada que deja fuera jun–sep |
|---|---|---|---|
| México | 70.5 % del total anual | 86.2 % | 18 % |
| CDMX | 75 % | 88 % | 15 % |
| Golfo, sureste y noreste | 52–63 % | 76–84 % | 25–33 % (lluvias de mayo, ciclones de oct–nov, nortes) |
| Pacífico occidental / monzón | 75–88 % | 79–97 % | 5–11 % |

Con 27 temporadas mensuales (1998–2024): el orden de los años casi no cambia (Spearman 0.70–0.94), pero la
clase seco/normal/húmedo cambia en 4–37 % de los años (nunca de seco a húmedo) y la significancia de la
tendencia cambia en Veracruz, Tabasco, Chiapas y BCS. **La ventana se fija antes de ver los resultados**:
elegirla según lo que "salga" infla el error tipo I.

**Por qué binario y no cuatro estaciones** (climatología mensual 1998–2024, % del total anual):

| | dic–feb | mar–may | jun–ago | sep–nov | lluvias | secas |
|---|---|---|---|---|---|---|
| México | 8 | 9 | 49 | 33 | 86 | 14 |
| CDMX | 3 | 13 | 55 | 29 | 88 | 12 |
| Tabasco | 14 | 9 | 35 | 42 | 78 | 22 |
| Baja California | **47** | 19 | 12 | 21 | 35 | **65** |

Los trimestres parten la temporada de lluvias a la mitad y dejan dos estaciones casi vacías: son un calendario
de temperatura de latitudes medias, no un régimen de lluvia. Quien los necesite (comparación internacional)
puede armarlos con los meses personalizados.

**Excepciones regionales**: Baja California tiene régimen de lluvias de **invierno** (65 % del total en la
temporada "seca" oficial); BCS reparte ~1/3 fuera de la oficial. Ahí "temporada de lluvias" no significa lo mismo.

**Datos mensuales** no se pueden cortar a medio mes; se usa la partición más cercana sin traslape:
- Lluvias → may–nov completo: **+5 %** de lluvia vs el corte exacto (de +1 a +11 % por estado). Sesgo casi constante:
  no afecta anomalías ni tendencias.
- Secas → dic–abr completo: **−34 %** (de −6 a −57 %): los días de transición (16–30 nov, 1–14 may) concentran mucha
  de la poca lluvia de secas. Para **totales** de secas usar medidas diarias (corte exacto).

**Comparaciones** (modo A | B): las cantidades (mm, días) se comparan en **% relativo**; las probabilidades y
porcentajes en **puntos porcentuales** (pp), porque un "% de un %" se malinterpreta fácilmente.

# Análisis estadístico — métodos y candados

`analysis.py` responde preguntas de precipitación con pruebas adecuadas a cómo es la lluvia
(muchos ceros, sesgo fuerte, estacionalidad y persistencia), y **se niega a reportar** lo que
los datos todavía no pueden sostener. Resultados en `data/analysis/`.

**Datos** (auditados: sin días ni meses faltantes):

| Escala | Producto | Periodo usado en estadística | Nota |
|---|---|---|---|
| Diario | GPM_3IMERGDF v07 (Final) | **2014-06-01 → 2025-09-30** | inicia con la era GPM; 2014 es año parcial |
| Diario | GPM_3IMERGDL v07 (Late) | — (sólo gráficas) | 2025-10 → hoy, preliminar sin calibrar con pluviómetros: **no entra a la estadística** |
| Mensual | GPM_3IMERGM v07 (Final) | **1998-01 → 2025-09** | oct–dic 2025 aún no publicados en Final |

El corte Final se lee de `is_late` (`analysis.final_end()`), así que avanza solo cuando se procese más Final.

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
| 31–57 % de días secos, asimetría 2.5–3.9 (2014–2025) | Nada de t-test ni OLS sobre mm. Pruebas de rangos (MK, Kruskal-Wallis), bootstrap, o gamma con masa en cero. |
| Totales mensuales: gamma gana por AIC | Respalda el **SPI** (gamma) para anomalías. Montos diarios: lognormal. |
| Estacionalidad (Kruskal-Wallis p < 0.01) | Comparar siempre contra el **mismo mes**; MK **estacional** para series mensuales. |
| Persistencia diaria r1 = 0.30–0.59 | Días no independientes: **bootstrap por bloques** (sino los IC salen demasiado estrechos → tipo I). |
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

Con 10 años (diario Final 2015–2024) sólo se detectarían cambios enormes (~75–110 %/década); con los 27 años
del mensual (1998–2024) se detectan ~18–27 %/década. **Las tendencias se prueban con el mensual.**

## 3. Preguntas y métricas mínimas viables

### ¿Cuándo viajar a Puerto Vallarta? (`travel`)
- **Métrica principal**: probabilidad de un *viaje seco* = ventana de N días (def. 7) con ≤ 1 día de lluvia.
  Secundarias: prob. diaria de lluvia (≥ 1 mm), de lluvia fuerte (≥ 10 mm), días de lluvia esperados.
- **Significancia**: bootstrap de la diferencia contra el mejor mes + Holm. Los meses *indistinguibles
  del mejor* forman la ventana recomendada; h de Cohen para ver si la diferencia importa.
- Resultado (Puerto Vallarta, diario Final 2014–2025, 11–12 años por mes → bootstrap **por años**, los IC ya
  incluyen variabilidad interanual): **abril** 100 % de viajes secos; **marzo (91 %) y mayo (94 %) indistinguibles**
  del mejor. Febrero (88 %) y nov–ene (~79–80 %) son estadísticamente peores (h ≈ 0.7–1.0). Jul–sep ≈ 0 %.
- Pendiente: sólo precipitación. Temperatura/humedad y temporada de huracanes (jun–nov) cambiarían el ranking.

### ¿Cuántas temporadas de lluvia en la ZMVM rebasaron el umbral de anomalía? (`season`)
- **Región**: CDMX + 17 municipios conurbados del Edomex, promedio ponderado por área
  (aproximación; la delimitación oficial tiene 76 municipios — editar `REGIONS["zmvm"]`).
- **Temporada**: oficial 15 may–15 nov por defecto (`--window lluvias|secas`; `--months` para otra). Ver §8.
- **Umbral**: SPI de la temporada con gamma ajustada a la referencia; rebase = SPI ≤ −1 (seco) o ≥ +1 (húmedo),
  esperado 15.9 % de los años cada uno. También Rx1day > P90 de la referencia.
- **Frecuencia**: k/n con IC exacto de Clopper-Pearson + prueba binomial contra lo esperado.
- Índices ETCCDI por temporada: total, días de lluvia, SDII, Rx1day, Rx5day, R20mm, R95p, CDD.
- Con el diario Final hay 11 temporadas (2014–2024): suficiente para los índices ETCCDI, **no** para una referencia
  de 20 años. Las anomalías y el umbral se responden con las medidas mensuales (`pct_normal`, `drought`),
  referencia 1998–2020.

### Tendencias (`trends`, `field`)
- Mensual: **MK estacional** (Hirsch 1982) y su versión correlacionada (Hirsch & Slack 1984).
- Anual y de temporada: **Hamed-Rao** (corrige autocorrelación) + MK original como contraste;
  si discrepan, manda Hamed-Rao. Pendiente de **Sen** (robusta a extremos).
- Campo: Hamed-Rao por celda + **FDR** → `field_trends.nc`, con el **mensual 1998–2024** (27 años).
  Con el diario (2015–2024, `--daily`) ~41 % de las celdas salen "significativas" y 99 % negativas: es la
  **variabilidad decenal** (2015–2018 húmedos, 2019–2023 secos; 2023 = 619 mm nacional), no una tendencia.
  Con 27 años: 10.7 % sin corregir y **520 celdas (3 %) con FDR**.

## 4. Regiones homogéneas y alto gradiente (`clusters`)
- K-means sobre la **climatología** de cada celda: forma del ciclo anual (12 fracciones mensuales),
  log del total anual, frecuencia de días con lluvia, SDII e índice de estacionalidad. k por silueta (hoy k = 6,
  con diario 2014–2025; separa un grupo de régimen invernal con pico en dic, el de Baja California).
- Se agrupa por **régimen**, no por tendencia: con pocos años, agrupar tendencias sería agrupar ruido.
  Después, la tendencia se prueba sobre la serie media de cada grupo (más potencia, menos pruebas).
- **Alto gradiente** = decil superior de |∇ log(total anual)| (montaña-costa, p. ej. Sierra Madre / Bahía de Banderas),
  y **frontera** = vecinos en otro grupo. Ahí la media del grupo no representa al sitio: usar celda o municipio.
- Silueta ~0.31: estructura moderada (los regímenes se mezclan gradualmente); tratar los bordes con cautela.

## 5. Qué falta para que las pruebas tengan potencia
1. **Mensual 1998-01 → 2025-09** (§6): referencia, SPI y tendencias con confianza alta.
   **Diario Final 2014-06 → 2025-09**: medidas diarias con confianza **media** (11–12 temporadas); pasan a alta
   con ≥ 15 (≈ 2029, o antes si se extiende el diario hacia la era TRMM, con cautela por el cambio de sensores).
2. Regenerar `zonal_series.nc` (`python zonal.py`) tras cada descarga del diario.
3. Validar IMERG contra estaciones del SMN en sitios de alto gradiente (IMERG subestima la lluvia orográfica).

---

## 6. Datos mensuales 1998 → 2025 (`download_monthly.py`)

**GPM_3IMERGM v07** (`C2723754851-GES_DISC`), misma malla de 0.1° que el diario. ~5 min de Harmony (1 job por año).

| Desbloquea (escala mensual) | No desbloquea (escala diaria) |
|---|---|
| Referencia 1998–2020 (23 años) para % de lo normal, terciles y SPI | Probabilidad de día con lluvia, rachas secas, Rx1day, "viaje seco" |
| SPI-1/3/6/12 (el SPI se define sobre totales mensuales) | Para estas: el **diario Final 2014-06 →** (confianza media) |
| Tendencias con 27 años (detectan ~20–27 %/década) | |
| Climatología robusta para `clusters` | |

- **Unidades**: el archivo trae `precipitation` en **mm/h** (tasa media); se convierte a mm/mes × 24 × días del mes.
- **Consistencia**: vs la suma de los diarios: r = 0.998 por celda-mes (2024); el total nacional mensual es
  **+2.4 a +4.6 %** mayor cada año 2014–2025, estable y sin deriva
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
| `wet_days` | ¿Probabilidad de que llueva un día? | diario | proporción ≥ umbral (1/10/20 mm), bootstrap por años | media |
| `dry_spell` | ¿Probabilidad de N días secos seguidos? | diario | ventanas sin lluvia, CDD por año | media |
| `intensity` | ¿Qué tan fuertes son los aguaceros? | diario | SDII, Rx1day, Rx5day | media |

Decisiones de diseño:
- **"% de aumento mes a mes" crudo no se ofrece**: el ciclo estacional lo domina (junio ≈ +200 % sobre mayo cada año).
  Se reporta el cambio observado **contra el cambio normal** y su percentil.
- **No significativo ≠ sin cambio**: `trend` reporta la tendencia mínima detectable para que la ausencia de
  evidencia no se lea como evidencia de ausencia (error tipo II).
- **Ocultar p-valores con confianza baja** evita que alguien encuentre "significancia" en 3 años de datos.
- Lugares: `{"cve": "09"}`, `{"point": [lat, lon]}` o `{"geometry": GeoJSON}` — los mismos que ya usa la app.

## 8. Temporadas: binario lluvias / secas, una ventana natural por escala

| Escala | Lluvias | Secas | Notas |
|---|---|---|---|
| **Diario** (2014-06→) | **15 may–15 nov** (oficial, exacto) | **16 nov–14 may** | corte al día |
| **Mensual** (1998→) | **may–nov** | **dic–abr** | meses completos |

En cada escala las dos temporadas son complementarias y **sin traslape**. Jun–sep queda como opción
"núcleo lluvioso" sólo para comparar con la literatura del monzón.

**Por qué no jun–sep** (diario Final 2015–2024, 10 años, `data/analysis/ventana_temporada.csv`):

| | Jun–sep | Oficial | Lluvia de temporada que deja fuera jun–sep |
|---|---|---|---|
| México | 68.0 % del total anual | 84.6 % | 20 % |
| CDMX | 72.9 % | 86.9 % | 16 % |
| Golfo, sureste y noreste | 50–63 % | 74–85 % | 24–36 % (lluvias de mayo, ciclones de oct–nov, nortes) |
| Pacífico occidental / monzón | 71–83 % | 76–93 % | 7–14 % |

Con 27 temporadas mensuales (1998–2024): el orden de los años casi no cambia (Spearman 0.70–0.94), pero la
clase seco/normal/húmedo cambia en 4–37 % de los años (nunca de seco a húmedo) y la significancia de la
tendencia cambia en Veracruz, Tabasco, Chiapas y BCS. **La ventana se fija antes de ver los resultados**:
elegirla según lo que "salga" infla el error tipo I.

**Por qué binario y no cuatro estaciones** (% del total anual; trimestres: mensual 1998–2024; lluvias/secas
oficiales: diario Final 2015–2024):

| | dic–feb | mar–may | jun–ago | sep–nov | lluvias | secas |
|---|---|---|---|---|---|---|
| México | 8 | 9 | 49 | 33 | 85 | 15 |
| CDMX | 3 | 13 | 55 | 29 | 87 | 13 |
| Tabasco | 14 | 9 | 35 | 42 | 80 | 20 |
| Baja California | **47** | 19 | 12 | 21 | 31 | **69** |

Los trimestres parten la temporada de lluvias a la mitad y dejan dos estaciones casi vacías: son un calendario
de temperatura de latitudes medias, no un régimen de lluvia. Quien los necesite (comparación internacional)
puede armarlos con los meses personalizados.

**Excepciones regionales**: Baja California tiene régimen de lluvias de **invierno** (69 % del total en la
temporada "seca" oficial); BCS reparte ~1/3 fuera de la oficial. Ahí "temporada de lluvias" no significa lo mismo.

**Contrapartida aceptada de los meses completos** (medida con el diario Final 2015–2024, 32 estados):
- Lluvias may–nov incluye las colas 1–14 may y 16–30 nov: **+4.5 %** de lluvia vs el oficial (de +1 a +19 %).
  Sesgo casi constante entre años: no afecta anomalías, % de lo normal, SPI ni tendencias.
- Secas dic–abr excluye esas colas: **−26 %** (de −9 a −45 %), porque concentran mucha de la poca lluvia
  de secas. Para **totales** de secas, preferir las medidas diarias.
- Por eso los totales de una escala y otra **no se comparan directamente**: la app nombra cada temporada
  con la ventana que usa ("15 may–15 nov" vs "may–nov").

**Comparaciones** (modo A | B): las cantidades (mm, días) se comparan en **% relativo**; las probabilidades y
porcentajes en **puntos porcentuales** (pp), porque un "% de un %" se malinterpreta fácilmente.

## 9. ¿Se puede integrar 2026 (IMERG Late, preliminar)? (`validate_late.py`)

Late (2025-10 → hoy) no tiene la calibración mensual con pluviómetros de Final. Se comparó Late vs Final
en dos años donde existen ambos (**oct-2023 → sep-2025**, 731 días, misma malla y máscara).
Resultados en `data/analysis/late_vs_final_*`.

**Alcance: compatible.** Misma malla 0.1°, unidades (mm/día), máscara y proceso; marcado día a día con `is_late`.
Cubre completa la temporada seca 2025–26 y el 99 % de la de lluvias 2026 (hasta 13-nov).

**Ocurrencia de lluvia (¿llovió o no?): compatible.**

| | Frecuencia ≥ 1 mm Final | Late | Dif. | POD | FAR | CSI |
|---|---|---|---|---|---|---|
| México | 21.2 % | 21.3 % | **+0.1 pp** | 0.88 | 0.12 | 0.79 |
| Regiones húmedas (grupos 5–6) | 29–33 % | 29–32 % | −1.0 a +0.4 pp | 0.89–0.92 | 0.08–0.10 | 0.83–0.84 |
| Regiones áridas (grupos 1–2) | 4–5 % | 5 % | +0.2 a +0.3 pp | 0.76 | 0.27–0.29 | 0.58–0.59 |

La frecuencia coincide dentro de ±1 pp en todas las regiones (margen ±5 pp). Las medidas de **ocurrencia**
(`wet_days`, `dry_spell`) podrían usar Late, marcado como preliminar; en zonas áridas con más cautela (CSI ~0.6).

**Cantidades (mm): NO compatibles.**
- Late es más húmedo: cociente regional Late/Final en lluvias **1.18** (IC95 1.07–1.24), 50 % central de las
  regiones-mes entre 0.97 y 1.45. Equivalencia ±10 % **rechazada** (TOST p = 0.83; secas p = 0.10).
  Nacional: Late > Final en +10 % o más en 14 de 24 meses; temporada de lluvias 2024 **+10.6 %**.
- Ese sesgo es **mayor que la variabilidad natural**: el CV interanual de la temporada de lluvias nacional es 9 %.
  Sumar Late haría que 2026 pareciera un año ~1 desviación estándar más húmedo sólo por el producto.
- Regionalmente es errático: temporadas de −47 % a +94 % según la región (las áridas, peor).
- **La corrección no transfiere**: factores Final/Late por celda × mes ajustados en el año 1 empeoran el año 2
  (error mediano 37 % → 43 %; sesgo nacional +13.5 % → +23.1 %). El sesgo cambia de un año a otro, así que no
  hay calibración confiable para 2026 con estos datos.

**Decisión (adoptada)**: IMERG Late (2025-10 → hoy) se usa **sólo en los mapas**, marcado como preliminar.
Series de tiempo, indicadores, medidas y análisis usan **sólo Final** (hasta `final_end()`, hoy 2025-09-30).
Aunque la ocurrencia de lluvia sí es equivalente, se excluye todo Late de la estadística para tener una sola
regla simple y evitar mezclar productos. Cuando NASA publique Final para esos meses (~3.5 meses de latencia),
`final_end()` avanza solo y esos meses entran a series y estadística sin cambios de código.

Nota de método: el TOST por **celda** (≈ 200 mil celdas-mes) se descartó: celdas vecinas comparten el clima del
día, así que no son independientes y cualquier sesgo medio pequeño "pasa". La unidad de prueba es región × mes.

# Fuente de los datos — Ejercicio práctico Tobit

## Conjunto de datos
**"Rain in Australia"** — observaciones meteorológicas diarias de estaciones de
Australia (~10 años, 2007–2017).

- **Repositorio público (Kaggle):** *Rain in Australia*, por Joe Young.
  https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package
- **Fuente original / autoridad de los datos:** **Australian Bureau of Meteorology
  (BoM)** — *Daily Weather Observations* y *Climate Data Online*.
  http://www.bom.gov.au/climate/data/  ·  http://www.bom.gov.au/climate/dwo/
- **Derechos:** © Commonwealth of Australia, Bureau of Meteorology.
  Datos publicados para uso público con atribución.

## Archivos
| Archivo | Descripción |
|---|---|
| `weatherAUS_raw.csv` | Datos crudos originales (todas las estaciones, 23 columnas, ~145.000 filas). |
| `lluvia_tobit.csv` | **Subconjunto curado** usado en el práctico (generado por `R/01_prep_lluvia.R`). |

## Subconjunto curado (`lluvia_tobit.csv`)
- **Estación:** Sydney (SydneyAirport).
- **Período:** 2009-01-01 a 2017-06-25.
- **Observaciones:** 2.948 días (filas con datos completos en las variables del modelo).
- **Censura:** 1.895 días con `lluvia_mm = 0` (**64,3 %**) → censura por la izquierda en 0,
  lo que motiva el uso del modelo **Tobit** en lugar de OLS.

### Variables
| Variable | Descripción | Unidad |
|---|---|---|
| `fecha` | Fecha de la observación | yyyy-mm-dd |
| `ciudad` | Estación meteorológica | — |
| `lluvia_mm` | **Variable dependiente:** precipitación del día (0 en días secos) | mm |
| `temp_max` | Temperatura máxima del día | °C |
| `viento_max` | Ráfaga de viento máxima | km/h |
| `humedad_pm` | Humedad relativa a las 15:00 | % |
| `presion_pm` | Presión atmosférica a las 15:00 | hPa |
| `nubosidad_pm` | Nubosidad a las 15:00 | oktas (0–8) |
# Fuente de los datos — Ejemplo B (modelo "bien comportado")

## Conjunto de datos
**Academic aptitude (`apt`)** — ejemplo didáctico clásico para enseñar el modelo Tobit,
publicado por el **UCLA Office of Advanced Research Computing (OARC / ex IDRE)**.

- **Tutorial desarrollado (R):** *Tobit Models — R Data Analysis Examples*, UCLA OARC.
  https://stats.oarc.ucla.edu/r/dae/tobit-models/
- **Archivo de datos (descarga directa):** https://stats.idre.ucla.edu/stat/data/tobit.csv
- **Naturaleza:** conjunto didáctico/ilustrativo ampliamente usado en la literatura para
  enseñar regresión censurada. Replica el comportamiento de una prueba estandarizada.

## Descripción
200 estudiantes rinden una prueba de aptitud académica (`apt`) cuyo puntaje va de 200 a
**800**. El puntaje **se topa en 800**: quienes "tocan techo" quedan registrados como 800
aunque su aptitud real sea mayor → **censura por la derecha** en 800 (17 estudiantes).

### Variables
| Variable | Descripción | Tipo |
|---|---|---|
| `id` | Identificador del estudiante | — |
| `apt` | **Variable dependiente:** puntaje de aptitud (censurado por arriba en 800) | 200–800 |
| `read` | Puntaje en lectura | continuo |
| `math` | Puntaje en matemática | continuo |
| `prog` | Tipo de programa: `academic`, `general`, `vocational` | categórica |

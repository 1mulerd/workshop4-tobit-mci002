# =====================================================================
#  Workshop 4 - MODELO TOBIT | EJEMPLO A (ingenieria): lluvia en Sydney
#  Caso: lluvia diaria (mm). La lluvia vale 0 en dias secos
#        => variable CENSURADA por la izquierda en 0 => se usa Tobit.
#  Es el caso "real y dificil": ajuste modesto y supuestos que NO se cumplen.
#
#  USO (RStudio):
#    1) Crear un proyecto en la carpeta raiz (donde estan data/ y R/).
#    2) Abrir este script y pulsar "Source" (o Ctrl+Shift+Enter).
# =====================================================================

# 0) Paquetes ----------------------------------------------------------
# Si no tienes los paquetes instalados, descomenta y corre esta linea primero:
# install.packages(c("censReg", "AER", "crch", "lmtest", "ggplot2"))

library(censReg)        # modelo Tobit por maxima verosimilitud
library(AER)            # tobit() alternativo (salida estilo libro)
library(crch)           # Tobit heterocedastico (test de homocedasticidad)
library(lmtest)         # test de razon de verosimilitud (lrtest)
library(ggplot2)        # graficos

set.seed(123)             # reproducibilidad

# 1) Cargar datos ------------------------------------------------------
#    Por defecto desde GitHub (link raw). Si tienes el repo clonado, usa
#    la copia local automaticamente.
url_github <- "https://raw.githubusercontent.com/1mulerd/workshop4-tobit-mci002/main/data/lluvia_tobit.csv"
ruta_local <- "data/lluvia_tobit.csv"
datos <- read.csv(if (file.exists(ruta_local)) ruta_local else url_github)

str(datos)

# 2) Analisis exploratorio (EDA) --------------------------------------
n        <- nrow(datos)
pct_cero <- mean(datos$lluvia_mm == 0)
cat(sprintf("\nObservaciones: %d | dias secos (lluvia=0): %d (%.1f%%)\n",
            n, sum(datos$lluvia_mm == 0), 100 * pct_cero))
summary(datos$lluvia_mm)

# Histograma: se ve la MASA de observaciones en 0 (lo que justifica Tobit)
g_hist <- ggplot(datos, aes(lluvia_mm)) +
  geom_histogram(binwidth = 2, fill = "#2c7fb8", color = "white") +
  labs(title = "Lluvia diaria en Sydney (mm)",
       subtitle = sprintf("%.0f%% de los dias = 0 mm (censura por la izquierda)", 100 * pct_cero),
       x = "Lluvia (mm)", y = "Frecuencia") +
  theme_minimal(base_size = 13)

# 3) Modelo INGENUO: regresion lineal (OLS) ----------------------------
#    Ignora la censura. Problema: puede predecir lluvias NEGATIVAS.
f <- lluvia_mm ~ temp_max + viento_max + humedad_pm + presion_pm + nubosidad_pm
m_ols <- lm(f, data = datos)

pred_ols <- fitted(m_ols)
cat(sprintf("\nOLS: %d predicciones NEGATIVAS de lluvia (imposibles fisicamente)\n",
            sum(pred_ols < 0)))

# 4) Modelo TOBIT (censurado por la izquierda en 0) --------------------
m_tobit <- censReg(f, left = 0, data = datos)   # estimacion por max. verosimilitud
summary(m_tobit)

# Cross-check con AER::tobit (mismos coeficientes, salida estilo "libro")
m_tobit_aer <- tobit(f, left = 0, data = datos)
summary(m_tobit_aer)

# Parametros: betas (efecto sobre la variable LATENTE y*) y sigma
b      <- coef(m_tobit)
sigma  <- exp(b["logSigma"])             # desviacion del error latente
beta   <- b[setdiff(names(b), "logSigma")]
cat(sprintf("\nSigma estimada (escala del error latente): %.3f\n", sigma))

# 5) Efectos marginales -------------------------------------------------
#    OJO: el coeficiente Tobit es el efecto sobre y* (latente), NO sobre la
#    lluvia observada. El efecto marginal sobre E[lluvia observada] es menor
#    (se "encoge" por la probabilidad de estar censurado).
cat("\n--- Efectos marginales sobre E[lluvia observada] ---\n")
print(summary(margEff(m_tobit)))

# 6) Metricas de ajuste y comparacion OLS vs Tobit ---------------------
# El R2 clasico NO es la vara correcta en Tobit. Reportamos 3 pseudo-R2:
X        <- model.matrix(f, data = datos)
xb       <- as.numeric(X %*% beta)              # indice latente predicho  x*beta
m_nulo   <- censReg(lluvia_mm ~ 1, left = 0, data = datos)

# (a) McKelvey-Zavoina: sobre la variable LATENTE y* (la mas apropiada en Tobit)
R2_mz    <- var(xb) / (var(xb) + sigma^2)
# (b) McFadden: basado en verosimilitud
R2_mcf   <- as.numeric(1 - logLik(m_tobit) / logLik(m_nulo))

# E[lluvia observada] segun Tobit:  E[y] = Phi(xb/s)*xb + s*phi(xb/s)  (siempre >= 0)
pred_tob <- pnorm(xb / sigma) * xb + sigma * dnorm(xb / sigma)
# (c) Correlacion^2 en la escala OBSERVADA (capacidad predictiva)
R2_obs   <- cor(datos$lluvia_mm, pred_tob)^2

cat(sprintf("\nPseudo-R2 -> McKelvey-Zavoina: %.3f (latente) | McFadden: %.3f | cor^2 obs.: %.3f\n",
            R2_mz, R2_mcf, R2_obs))

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
r2obs <- function(obs, pred) cor(obs, pred)^2   # R2 en la escala observada (capac. predictiva)
tabla <- data.frame(
  Modelo          = c("OLS (lm)", "Tobit"),
  LogLik          = c(as.numeric(logLik(m_ols)),   as.numeric(logLik(m_tobit))),
  AIC             = c(AIC(m_ols),                  AIC(m_tobit)),
  RMSE            = c(rmse(datos$lluvia_mm, pred_ols), rmse(datos$lluvia_mm, pred_tob)),
  R2_obs          = c(r2obs(datos$lluvia_mm, pred_ols), r2obs(datos$lluvia_mm, pred_tob)),
  Pred_negativas  = c(sum(pred_ols < 0),           sum(pred_tob < 0))
)
cat("\n--- Comparacion de modelos ---\n"); print(tabla, row.names = FALSE, digits = 4)

# 7) DIAGNOSTICO DE SUPUESTOS (clave en Tobit) -------------------------
#    A diferencia de OLS, si en Tobit se violan la NORMALIDAD o la
#    HOMOCEDASTICIDAD del error latente, el estimador MLE es SESGADO e
#    INCONSISTENTE (no solo ineficiente). Por eso hay que testearlas.

# --- 7.1 Residuos generalizados (Chesher-Irish) -----------------------
#   No censurado (y>0):  (y - xb)/sigma
#   Censurado  (y=0):   -phi(xb/sigma) / Phi(-xb/sigma)   (E[error | y*<=0])
cens   <- datos$lluvia_mm == 0
res_g  <- ifelse(cens,
                 -dnorm(xb / sigma) / pnorm(-xb / sigma),
                 (datos$lluvia_mm - xb) / sigma)

# --- 7.2 Normalidad: asimetria, curtosis y Jarque-Bera ----------------
n   <- length(res_g)
S   <- mean((res_g - mean(res_g))^3) / sd(res_g)^3          # asimetria
K   <- mean((res_g - mean(res_g))^4) / sd(res_g)^4          # curtosis
JB  <- n / 6 * (S^2 + (K - 3)^2 / 4)                        # ~ chi2(2) bajo normalidad
cat(sprintf("\n[Normalidad] Asimetria=%.2f  Curtosis=%.2f  Jarque-Bera=%.1f  p=%.3g\n",
            S, K, JB, pchisq(JB, df = 2, lower.tail = FALSE)))

# Evidencia adicional: comparar distribuciones del error (AIC). Si la t de
# Student (colas pesadas) ajusta mucho mejor, el supuesto normal es dudoso.
m_norm <- crch(f, data = datos, left = 0, dist = "gaussian")
m_t    <- crch(f, data = datos, left = 0, dist = "student")
cat(sprintf("[Normalidad] AIC error normal=%.0f  vs  AIC error t-Student=%.0f  (menor=mejor)\n",
            AIC(m_norm), AIC(m_t)))

# --- 7.3 Homocedasticidad: test LR (Tobit homo vs heterocedastico) ----
#   m_homo: escala (log sigma) constante.  m_het: escala ~ predictores.
m_homo <- crch(f, data = datos, left = 0, dist = "gaussian")
m_het  <- crch(lluvia_mm ~ temp_max + viento_max + humedad_pm + presion_pm + nubosidad_pm |
                           temp_max + viento_max + humedad_pm + presion_pm + nubosidad_pm,
               data = datos, left = 0, dist = "gaussian")
cat("\n[Homocedasticidad] Test de razon de verosimilitud (H0: errores homocedasticos):\n")
print(lrtest(m_homo, m_het))

# --- 7.4 QQ-plot de los residuos generalizados ------------------------
g_qq <- ggplot(data.frame(r = res_g), aes(sample = r)) +
  stat_qq(alpha = 0.3, color = "#2c7fb8") + stat_qq_line(color = "red", linetype = 2) +
  labs(title = "QQ-plot de residuos generalizados",
       subtitle = "Desvios de la recta = alejamiento de la normalidad",
       x = "Cuantiles teoricos N(0,1)", y = "Cuantiles muestrales") +
  theme_minimal(base_size = 13)

# 8) Visualizaciones ----------------------------------------------------
# (a) Predicho (Tobit) vs observado
g_fit <- ggplot(data.frame(obs = datos$lluvia_mm, pred = pred_tob),
                aes(pred, obs)) +
  geom_point(alpha = 0.25, color = "#2c7fb8") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
  labs(title = "Tobit: lluvia predicha vs observada",
       x = "E[lluvia] predicha (mm)", y = "Lluvia observada (mm)") +
  theme_minimal(base_size = 13)

# (b) Efecto parcial de la humedad sobre E[lluvia] (resto en su media)
rejilla <- data.frame(
  temp_max     = mean(datos$temp_max),
  viento_max   = mean(datos$viento_max),
  humedad_pm   = seq(min(datos$humedad_pm), max(datos$humedad_pm), length.out = 100),
  presion_pm   = mean(datos$presion_pm),
  nubosidad_pm = mean(datos$nubosidad_pm)
)
xb_g <- as.numeric(model.matrix(delete.response(terms(f)), rejilla) %*% beta)
rejilla$Ey <- pnorm(xb_g / sigma) * xb_g + sigma * dnorm(xb_g / sigma)

g_efecto <- ggplot(rejilla, aes(humedad_pm, Ey)) +
  geom_line(color = "#d95f0e", linewidth = 1.2) +
  labs(title = "Efecto parcial de la humedad sobre la lluvia esperada",
       subtitle = "Resto de variables fijadas en su media",
       x = "Humedad a las 15 h (%)", y = "E[lluvia] (mm)") +
  theme_minimal(base_size = 13)

# 9) Guardar figuras ----------------------------
if (!dir.exists("img2")) dir.create("img2", recursive = TRUE)
ggsave("img2/lluvia_hist.png",        g_hist,   width = 7, height = 4.2, dpi = 150)
ggsave("img2/lluvia_pred_vs_obs.png", g_fit,    width = 7, height = 4.2, dpi = 150)
ggsave("img2/lluvia_efecto.png",      g_efecto, width = 7, height = 4.2, dpi = 150)
ggsave("img2/lluvia_qq.png",          g_qq,     width = 7, height = 4.2, dpi = 150)
cat("\nFiguras guardadas en img2/. Listo.\n")

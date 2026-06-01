# =====================================================================
#  Workshop 4 - MODELO TOBIT | EJEMPLO B: aptitud academica (apt)
#  Caso: puntaje de aptitud (200-800) que se TOPA en 800.
#        => variable CENSURADA por la DERECHA en 800 => se usa Tobit.
#  Es el caso "de manual": buen ajuste y supuestos que SI se cumplen.
#  (Contrapunto del ejemplo de lluvia.)
#
#  USO (RStudio):
#    1) Crear un proyecto en la carpeta raiz (donde estan data/ y R/).
#    2) Abrir este script y pulsar "Source" (o Ctrl+Shift+Enter).
# =====================================================================

# 0) Paquetes ----------------------------------------------------------
# Si no tienes los paquetes instalados, descomenta y corre esta linea primero:
# install.packages(c("censReg", "AER", "crch", "lmtest", "ggplot2"))

library(censReg)        # Tobit por maxima verosimilitud
library(AER)            # tobit() alternativo (salida estilo libro)
library(crch)           # Tobit heterocedastico (test de homocedasticidad)
library(lmtest)         # test de razon de verosimilitud (lrtest)
library(ggplot2)        # graficos
set.seed(123)

U <- 800   # limite de censura por la derecha (techo del puntaje)

# 1) Cargar datos ------------------------------------------------------
#    Por defecto desde GitHub (raw); usa copia local si el repo esta clonado.
url_github <- "https://raw.githubusercontent.com/1mulerd/workshop4-tobit-mci002/main/data/apt.csv"
ruta_local <- "data/apt.csv"
datos <- read.csv(if (file.exists(ruta_local)) ruta_local else url_github)
datos$prog <- factor(datos$prog)          # tipo de programa = categorica
str(datos)

# 2) Analisis exploratorio (EDA) --------------------------------------
pct_top <- mean(datos$apt == U)
cat(sprintf("\nObservaciones: %d | en el techo (apt=%d): %d (%.1f%%)\n",
            nrow(datos), U, sum(datos$apt == U), 100 * pct_top))

g_hist <- ggplot(datos, aes(apt)) +
  geom_histogram(binwidth = 20, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = U, color = "#d95f0e", linetype = 2, linewidth = 1) +
  annotate("text", x = U, y = Inf, label = "techo = 800", vjust = 1.5, hjust = 1.1,
           color = "#d95f0e") +
  labs(title = "Puntaje de aptitud (apt)",
       subtitle = sprintf("%.0f%% de los alumnos topan en 800 (censura por la derecha)", 100 * pct_top),
       x = "Aptitud (apt)", y = "Frecuencia") +
  theme_minimal(base_size = 13)

# 3) Modelo INGENUO: OLS ----------------------------------------------
#    Ignora la censura. Problema: puede predecir puntajes IMPOSIBLES (> 800).
f <- apt ~ read + math + prog
m_ols    <- lm(f, data = datos)
pred_ols <- fitted(m_ols)
cat(sprintf("\nOLS: %d predicciones POR ENCIMA del techo (apt > %d, imposibles)\n",
            sum(pred_ols > U), U))

# 4) Modelo TOBIT (censurado por la derecha en 800) --------------------
m_tobit <- censReg(f, left = -Inf, right = U, data = datos)   # max. verosimilitud
summary(m_tobit)
m_tobit_aer <- tobit(f, left = -Inf, right = U, data = datos) # cross-check
summary(m_tobit_aer)

b     <- coef(m_tobit)
sigma <- exp(b["logSigma"])
beta  <- b[setdiff(names(b), "logSigma")]
cat(sprintf("\nSigma estimada (escala del error latente): %.2f\n", sigma))

# 5) Efectos marginales sobre E[apt observado] -------------------------
cat("\n--- Efectos marginales sobre E[apt observado] ---\n")
print(summary(margEff(m_tobit)))

# 6) Metricas de ajuste y comparacion OLS vs Tobit ---------------------
X      <- model.matrix(f, data = datos)
xb     <- as.numeric(X %*% beta)                 # indice latente predicho
m_nulo <- censReg(apt ~ 1, left = -Inf, right = U, data = datos)

R2_mz  <- var(xb) / (var(xb) + sigma^2)          # McKelvey-Zavoina (latente) <- principal
R2_mcf <- as.numeric(1 - logLik(m_tobit) / logLik(m_nulo))   # McFadden

# E[apt observado] con censura por la derecha en U:
#   z = (U - xb)/s ;  E[y] = xb*Phi(z) - s*phi(z) + U*(1 - Phi(z))   (siempre <= U)
z        <- (U - xb) / sigma
pred_tob <- xb * pnorm(z) - sigma * dnorm(z) + U * (1 - pnorm(z))
R2_obs   <- cor(datos$apt, pred_tob)^2

cat(sprintf("\nPseudo-R2 -> McKelvey-Zavoina: %.3f (latente) | McFadden: %.3f | cor^2 obs.: %.3f\n",
            R2_mz, R2_mcf, R2_obs))

rmse  <- function(obs, pred) sqrt(mean((obs - pred)^2))
r2obs <- function(obs, pred) cor(obs, pred)^2
tabla <- data.frame(
  Modelo         = c("OLS (lm)", "Tobit"),
  LogLik         = c(as.numeric(logLik(m_ols)),  as.numeric(logLik(m_tobit))),
  AIC            = c(AIC(m_ols),                 AIC(m_tobit)),
  RMSE           = c(rmse(datos$apt, pred_ols),  rmse(datos$apt, pred_tob)),
  R2_obs         = c(r2obs(datos$apt, pred_ols), r2obs(datos$apt, pred_tob)),
  Pred_sobre_800 = c(sum(pred_ols > U),          sum(pred_tob > U))
)
cat("\n--- Comparacion de modelos ---\n"); print(tabla, row.names = FALSE, digits = 4)

# 7) DIAGNOSTICO DE SUPUESTOS ------------------------------------------
#    En Tobit, violar normalidad/homocedasticidad => estimador SESGADO e
#    INCONSISTENTE. Aqui esperamos que SI se cumplan (a diferencia de la lluvia).

# 7.1 Residuos generalizados (censura por la derecha)
#    No censurado (y<U): (y - xb)/sigma
#    Censurado  (y=U):    phi(z)/(1 - Phi(z))   con z=(U-xb)/sigma
cens  <- datos$apt == U
res_g <- ifelse(cens, dnorm(z) / (1 - pnorm(z)), (datos$apt - xb) / sigma)

# 7.2 Normalidad: asimetria, curtosis y Jarque-Bera
n  <- length(res_g)
S  <- mean((res_g - mean(res_g))^3) / sd(res_g)^3
K  <- mean((res_g - mean(res_g))^4) / sd(res_g)^4
JB <- n / 6 * (S^2 + (K - 3)^2 / 4)
cat(sprintf("\n[Normalidad] Asimetria=%.2f  Curtosis=%.2f  Jarque-Bera=%.1f  p=%.3g\n",
            S, K, JB, pchisq(JB, df = 2, lower.tail = FALSE)))
m_norm <- crch(f, data = datos, left = -Inf, right = U, dist = "gaussian")
m_t    <- crch(f, data = datos, left = -Inf, right = U, dist = "student")
cat(sprintf("[Normalidad] AIC error normal=%.0f  vs  AIC error t-Student=%.0f  (menor=mejor)\n",
            AIC(m_norm), AIC(m_t)))

# 7.3 Homocedasticidad: test LR (Tobit homo vs heterocedastico)
m_homo <- crch(f, data = datos, left = -Inf, right = U, dist = "gaussian")
m_het  <- crch(apt ~ read + math + prog | read + math + prog,
               data = datos, left = -Inf, right = U, dist = "gaussian")
cat("\n[Homocedasticidad] Test de razon de verosimilitud (H0: errores homocedasticos):\n")
print(lrtest(m_homo, m_het))

# 7.4 QQ-plot de residuos generalizados
g_qq <- ggplot(data.frame(r = res_g), aes(sample = r)) +
  stat_qq(alpha = 0.4, color = "#2c7fb8") + stat_qq_line(color = "red", linetype = 2) +
  labs(title = "QQ-plot de residuos generalizados (apt)",
       subtitle = "Puntos sobre la recta = normalidad razonable",
       x = "Cuantiles teoricos N(0,1)", y = "Cuantiles muestrales") +
  theme_minimal(base_size = 13)

# 8) Visualizaciones ----------------------------------------------------
g_fit <- ggplot(data.frame(obs = datos$apt, pred = pred_tob), aes(pred, obs)) +
  geom_point(alpha = 0.5, color = "#2c7fb8") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
  labs(title = "Tobit: aptitud predicha vs observada",
       x = "E[apt] predicha", y = "apt observada") +
  theme_minimal(base_size = 13)

rejilla <- data.frame(read = mean(datos$read),
                      math = seq(min(datos$math), max(datos$math), length.out = 100),
                      prog = factor("academic", levels = levels(datos$prog)))
xb_g <- as.numeric(model.matrix(delete.response(terms(f)), rejilla) %*% beta)
zg   <- (U - xb_g) / sigma
rejilla$Ey <- xb_g * pnorm(zg) - sigma * dnorm(zg) + U * (1 - pnorm(zg))
g_efecto <- ggplot(rejilla, aes(math, Ey)) +
  geom_line(color = "#d95f0e", linewidth = 1.2) +
  labs(title = "Efecto parcial de 'math' sobre la aptitud esperada",
       subtitle = "read en su media; prog = academic",
       x = "Puntaje en matematica", y = "E[apt]") +
  theme_minimal(base_size = 13)

# 9) Guardar figuras ----------------------------------------------------
if (!dir.exists("img3")) dir.create("img3", recursive = TRUE)
ggsave("img3/apt_hist.png",        g_hist,   width = 7, height = 4.2, dpi = 150)
ggsave("img3/apt_pred_vs_obs.png", g_fit,    width = 7, height = 4.2, dpi = 150)
ggsave("img3/apt_efecto.png",      g_efecto, width = 7, height = 4.2, dpi = 150)
ggsave("img3/apt_qq.png",          g_qq,     width = 7, height = 4.2, dpi = 150)
cat("\nFiguras guardadas en img3/. Listo.\n")

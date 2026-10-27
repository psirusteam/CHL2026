##############################################################################
## Dia 1 - Sesion 2: Fundamentos de la inferencia bayesiana I              ##
## Modelos uniparametricos: Bernoulli-Beta, Poisson-Gamma y Normal.        ##
## Ejecutar con el directorio de trabajo en Dia1_S2_Bayes_uniparametrico/: ##
##   source("Rcode/Bayes_uniparametrico.R")                                 ##
##############################################################################

source("Rcode/utilidades.R")

## ---- paquetes ----
library(tidyverse)
library(rstan)
library(bayesplot)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_minimal(base_size = 15))

## ---- datos ----
segmentos <- readRDS("data/segmentos_aysen.rds")

tabla_datos <- segmentos %>%
  count(zona, clasificacion_cobertura) %>%
  pivot_wider(names_from = clasificacion_cobertura, values_from = n)
tabla_datos

## ---- guardar-datos ----
guardar(tabla_datos, "tabla_datos")
guardar_valor("n_segmentos", nrow(segmentos))

# =============================================================================
# Modelo Bernoulli-Beta
# =============================================================================

## ---- bernoulli-analitica ----
y <- as.integer(segmentos$clasificacion_cobertura == "Completa")
n <- length(y)
a <- 1; b <- 1                          # previa uniforme Beta(1, 1)

a_post <- a + sum(y)
b_post <- b + n - sum(y)

bern_analitica <- tibble(
  n = n, completos = sum(y),
  media_post = a_post / (a_post + b_post),
  desv_post = sqrt(a_post * b_post / ((a_post + b_post)^2 * (a_post + b_post + 1))),
  ic_025 = qbeta(0.025, a_post, b_post),
  ic_975 = qbeta(0.975, a_post, b_post)
)
bern_analitica

## ---- guardar-bernoulli-analitica ----
guardar(bern_analitica, "bern_analitica")
guardar_valor("bern_completos", sum(y))

## ---- bernoulli-previas ----
set.seed(123)
muestra <- sample(y, 20)
previas <- tibble(previa = c("Beta(1, 1): uniforme", "Beta(2, 8): pesimista", "Beta(8, 2): optimista"),
                  a = c(1, 2, 8), b = c(1, 8, 2))
grilla <- seq(0.001, 0.999, length.out = 400)

curvas <- crossing(previas, datos = c("20 segmentos", paste(n, "segmentos"))) %>%
  mutate(y_obs = if_else(datos == "20 segmentos", sum(muestra), sum(y)),
         n_obs = if_else(datos == "20 segmentos", 20L, n)) %>%
  rowwise() %>%
  reframe(previa, datos, theta = grilla,
          Previa = dbeta(grilla, a, b),
          Posterior = dbeta(grilla, a + y_obs, b + n_obs - y_obs)) %>%
  pivot_longer(c(Previa, Posterior), names_to = "tipo", values_to = "densidad") %>%
  mutate(datos = factor(datos, levels = c("20 segmentos", paste(n, "segmentos"))))

grafico_previas <- ggplot(curvas, aes(theta, densidad, colour = previa, linetype = tipo)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ datos, scales = "free_y") +
  labs(x = expression(theta), y = "Densidad", colour = NULL, linetype = NULL) +
  theme(legend.position = "bottom")
guardar_figura(grafico_previas, "grafico_previas")

tabla_previas <- previas %>%
  crossing(tibble(datos = c("20 segmentos", "todos"), y_obs = c(sum(muestra), sum(y)), n_obs = c(20, n))) %>%
  mutate(media_posterior = (a + y_obs) / (a + b + n_obs),
         peso_datos = n_obs / (a + b + n_obs)) %>%
  select(previa, datos, media_posterior, peso_datos)
guardar(tabla_previas, "tabla_previas")

## ---- bernoulli-stan ----
modelo_bern <- stan_model("stan/bernoulli_beta.stan")

datos_bern <- list(n = n, y = y, a = 1, b = 1)

fit_bern <- sampling(
  modelo_bern, data = datos_bern,
  chains = 4, iter = 2000, warmup = 1000,
  seed = 1234, refresh = 0
)

## ---- bernoulli-resumen ----
resumen_bern <- resumen_stan(fit_bern, "theta")
resumen_bern

## ---- guardar-bernoulli-stan ----
guardar(resumen_bern, "resumen_bern")

## ---- bernoulli-comparacion ----
theta_draws <- as.data.frame(fit_bern)$theta
grafico_bern <- ggplot(tibble(theta = theta_draws), aes(theta)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#9ecae1", colour = "white") +
  stat_function(fun = dbeta, args = list(shape1 = a_post, shape2 = b_post),
                colour = "#c00000", linewidth = 1.2) +
  labs(x = expression(theta), y = "Densidad",
       subtitle = "Histograma: draws de Stan · Línea roja: posterior Beta exacta")
guardar_figura(grafico_bern, "grafico_bern", alto = 4.8)

## ---- bernoulli-predictiva ----
# Predictiva: numero de segmentos completos entre m = 20 segmentos nuevos.
# Exacta: Beta-Binomial(m, a_post, b_post). Simulada: Binomial(m, theta_s).
m <- 20
dbetabinom <- function(k, m, a, b) exp(lchoose(m, k) + lbeta(k + a, m - k + b) - lbeta(a, b))
predictiva <- tibble(k = 0:m,
                     exacta = dbetabinom(0:m, m, a_post, b_post),
                     simulada = as.numeric(table(factor(rbinom(length(theta_draws), m, theta_draws),
                                                        levels = 0:m))) / length(theta_draws))
grafico_predictiva <- predictiva %>%
  ggplot(aes(k)) +
  geom_col(aes(y = simulada), fill = "#9ecae1") +
  geom_point(aes(y = exacta), colour = "#c00000", size = 2.5) +
  labs(x = "Segmentos completos entre 20 nuevos", y = "Probabilidad",
       subtitle = "Barras: simulación con los draws · Puntos: Beta-Binomial exacta")
guardar_figura(grafico_predictiva, "grafico_predictiva_bern", alto = 4.5)

# =============================================================================
# Modelo Poisson-Gamma con exposicion
# =============================================================================

## ---- poisson-datos ----
completos <- segmentos %>%
  filter(clasificacion_cobertura == "Completa", !excluido_colectivo, n_vpo > 0)

a <- 0.01; b <- 0.01                  # previa vaga: media 1, varianza 100
a_post <- a + sum(completos$n_per)
b_post <- b + sum(completos$n_vpo)

pg_analitica <- tibble(
  segmentos = nrow(completos),
  personas = sum(completos$n_per), viviendas = sum(completos$n_vpo),
  media_post = a_post / b_post, desv_post = sqrt(a_post) / b_post,
  ic_025 = qgamma(0.025, a_post, b_post),
  ic_975 = qgamma(0.975, a_post, b_post)
)
pg_analitica

## ---- guardar-poisson-analitica ----
guardar(pg_analitica, "pg_analitica")

## ---- poisson-sensibilidad ----
sensibilidad <- tibble(previa = c("Gamma(0,01; 0,01): vaga", "Gamma(30; 10): media 3",
                                  "Gamma(300; 100): media 3, muy informativa", "Gamma(2; 1): media 2"),
                       a = c(0.01, 30, 300, 2), b = c(0.01, 10, 100, 1)) %>%
  mutate(media_previa = a / b,
         media_posterior = (a + sum(completos$n_per)) / (b + sum(completos$n_vpo)),
         peso_datos = sum(completos$n_vpo) / (b + sum(completos$n_vpo)))
guardar(select(sensibilidad, -a, -b), "tabla_sensibilidad")

## ---- poisson-stan ----
modelo_pg <- stan_model("stan/poisson_gamma.stan")

datos_pg <- list(
  n = nrow(completos),
  y = completos$n_per,
  v = completos$n_vpo,
  a = 0.01, b = 0.01
)

fit_pg <- sampling(
  modelo_pg, data = datos_pg,
  chains = 4, iter = 2000, warmup = 1000,
  seed = 1234, refresh = 0
)

## ---- poisson-resumen ----
resumen_pg <- resumen_stan(fit_pg, c("delta", "y_nuevo"))
resumen_pg

## ---- guardar-poisson-stan ----
guardar(resumen_pg, "resumen_pg")

## ---- poisson-prediccion ----
y_nuevo <- as.data.frame(fit_pg)$y_nuevo

# Predictiva exacta: Binomial Negativa con r = a_post y p = b_post / (b_post + 50)
prediccion <- tibble(y = 60:200) %>%
  mutate(exacta = dnbinom(y, size = a_post, prob = b_post / (b_post + 50)))

grafico_prediccion <- ggplot() +
  geom_histogram(aes(y_nuevo, after_stat(density)), binwidth = 2, fill = "#2e75b6", colour = "white") +
  geom_line(aes(y, exacta), data = prediccion, colour = "#c00000", linewidth = 1.1) +
  coord_cartesian(xlim = c(80, 175)) +
  labs(x = "Personas en un segmento nuevo de 50 viviendas", y = "Densidad",
       subtitle = "Histograma: generated quantities · Línea roja: Binomial Negativa exacta")
guardar_figura(grafico_prediccion, "grafico_prediccion", alto = 4.5)
guardar(tibble(cuantil = c("2,5%", "50%", "97,5%"), y_nuevo = quantile(y_nuevo, c(0.025, 0.5, 0.975))),
        "cuantiles_y_nuevo")

## ---- poisson-ppc ----
delta_draws <- as.data.frame(fit_pg)$delta
dispersion <- function(y, mu) sum((y - mu)^2 / mu) / (length(y) - 1)

set.seed(1)
phi_rep <- map_dbl(sample(delta_draws, 500), function(d) {
  mu <- d * completos$n_vpo
  dispersion(rpois(length(mu), mu), mu)
})
phi_obs <- dispersion(completos$n_per, mean(delta_draws) * completos$n_vpo)

c(observado = phi_obs, replicas_q95 = quantile(phi_rep, 0.95),
  p_bayes = mean(phi_rep >= phi_obs))

## ---- guardar-poisson-ppc ----
guardar_valor("phi_obs", phi_obs)
guardar_valor("phi_rep_q95", unname(quantile(phi_rep, 0.95)))
grafico_ppc <- ggplot(tibble(phi_rep), aes(phi_rep)) +
  geom_histogram(bins = 40, fill = "#9ecae1", colour = "white") +
  geom_vline(xintercept = phi_obs, colour = "#c00000", linewidth = 1.2) +
  labs(x = "Índice de dispersión de Pearson", y = "Réplicas",
       subtitle = "Histograma: réplicas del modelo · Línea roja: datos observados")
guardar_figura(grafico_ppc, "grafico_ppc", alto = 4.2)

## ---- poisson-tamano ----
d_hat <- mean(delta_draws)
tabla_tamano <- completos %>%
  mutate(viviendas = cut(n_vpo, c(0, 10, 30, 60, Inf),
                         labels = c("1–10", "11–30", "31–60", "más de 60"))) %>%
  group_by(viviendas) %>%
  summarise(segmentos = n(),
            de_observada = sd(n_per / n_vpo),
            de_poisson = mean(sqrt(d_hat / n_vpo))) %>%
  mutate(razon = de_observada / de_poisson)
tabla_tamano

## ---- guardar-poisson-tamano ----
guardar(tabla_tamano, "tabla_tamano")

# =============================================================================
# Modelo Normal con varianza conocida
# =============================================================================

## ---- normal-analitica ----
aptos <- segmentos %>% filter(apto_modelo)
y_log <- log(aptos$delta)
sigma <- sd(y_log)                        # aproximacion didactica: se fija, no es conocida externamente
mu0 <- log(3); tau0 <- 0.5                # previa: alrededor de 3 personas por vivienda

prec_n <- 1 / tau0^2 + length(y_log) / sigma^2
mu_n <- (mu0 / tau0^2 + sum(y_log) / sigma^2) / prec_n
normal_analitica <- tibble(media_post = mu_n, desv_post = sqrt(1 / prec_n),
                           media_muestral = mean(y_log), n = length(y_log), sigma = sigma)
normal_analitica

## ---- guardar-normal-analitica ----
guardar(normal_analitica, "normal_analitica")

## ---- normal-stan ----
modelo_normal <- stan_model("stan/normal_media.stan")
fit_normal <- sampling(
  modelo_normal,
  data = list(n = length(y_log), y = y_log, sigma = sigma, mu0 = mu0, tau0 = tau0),
  chains = 4, iter = 2000, seed = 1234, refresh = 0
)
resumen_normal <- resumen_stan(fit_normal, "mu")

## ---- guardar-normal-stan ----
guardar(resumen_normal, "resumen_normal")

## ---- normal-contraccion ----
set.seed(2026)
tamanos <- c(2, 5, 10, 25, 50, 100, 250, 500, length(y_log))
contraccion <- map_dfr(tamanos, function(m) {
  ym <- sample(y_log, m)
  prec <- 1 / tau0^2 + m / sigma^2
  tibble(n = m,
         media = (mu0 / tau0^2 + sum(ym) / sigma^2) / prec,
         li = media - 1.96 / sqrt(prec), ls = media + 1.96 / sqrt(prec),
         peso_datos = (m / sigma^2) / prec)
})

grafico_contraccion <- ggplot(contraccion, aes(n, exp(media))) +
  geom_hline(yintercept = exp(mean(y_log)), linetype = "dashed", colour = "#c00000", linewidth = 1) +
  geom_hline(yintercept = exp(mu0), linetype = "dotted") +
  geom_pointrange(aes(ymin = exp(li), ymax = exp(ls)), colour = "#1f4e79") +
  geom_text(aes(label = paste0(round(100 * peso_datos), "%")), vjust = -1.2, size = 4.5) +
  scale_x_log10(breaks = tamanos) +
  labs(x = "Segmentos usados (escala log)", y = "exp(μ): personas por vivienda",
       subtitle = "Etiqueta: peso de los datos · roja: exp(media de log-densidad) · punteada: exp(media previa)")
guardar_figura(grafico_contraccion, "grafico_contraccion", alto = 5)


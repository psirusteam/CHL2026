##############################################################################
## Dia 1 - Sesion 3: Fundamentos de la inferencia bayesiana II             ##
## Modelos multiparametricos y metodos MCMC.                                ##
## Ejecutar con el directorio de trabajo en Dia1_S3_Bayes_multiparametrico/##
##   source("Rcode/Bayes_multiparametrico.R")                               ##
##############################################################################

source("Rcode/utilidades.R")

## ---- paquetes ----
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)
library(patchwork)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_minimal(base_size = 15))
color_scheme_set("blue")

## ---- datos ----
segmentos <- readRDS("data/segmentos_aysen.rds")
aptos <- segmentos %>% filter(apto_modelo)
stopifnot(nrow(aptos) > 2, all(is.finite(aptos$delta)), all(aptos$delta > 0))
y <- log(aptos$delta)
stopifnot(is.finite(var(y)), var(y) > 0)

# =============================================================================
# Normal con media y varianza desconocidas
# =============================================================================

## ---- distribucion-delta ----
grafico_log <- bind_rows(tibble(escala = "δ: personas por vivienda", valor = aptos$delta),
                         tibble(escala = "log δ", valor = log(aptos$delta))) %>%
  ggplot(aes(valor)) +
  geom_histogram(bins = 40, fill = "#2e75b6", colour = "white") +
  facet_wrap(~ escala, scales = "free") +
  labs(x = NULL, y = "Segmentos")
guardar_figura(grafico_log, "grafico_log", alto = 4.5)
asimetria <- function(x) mean((x - mean(x))^3) / sd(x)^3
guardar(tibble(escala = c("delta", "log(delta)"),
               asimetria = c(asimetria(aptos$delta), asimetria(log(aptos$delta)))),
        "tabla_asimetria")

## ---- normal-analitica ----
# Con la previa no informativa p(mu, sigma^2) proporcional a 1/sigma^2:
#   mu | y ~ t_{n-1}(ybar, s^2/n)     sigma^2 | y ~ Inv-chi^2(n-1, s^2)
n <- length(y); ybar <- mean(y); s2 <- var(y)
normal_exacta <- tibble(
  parametro = c("mu", "sigma"),
  media = c(ybar, sqrt(s2) * sqrt((n - 1) / 2) * exp(lgamma((n - 2) / 2) - lgamma((n - 1) / 2))),
  q2.5 = c(ybar + qt(0.025, n - 1) * sqrt(s2 / n), sqrt((n - 1) * s2 / qchisq(0.975, n - 1))),
  q97.5 = c(ybar + qt(0.975, n - 1) * sqrt(s2 / n), sqrt((n - 1) * s2 / qchisq(0.025, n - 1)))
)
normal_exacta

## ---- guardar-normal-analitica ----
guardar(normal_exacta, "normal_exacta")
guardar_valor("n_aptos", n)

## ---- normal-stan ----
modelo_normal <- stan_model("stan/normal_media_varianza.stan")

datos_normal <- list(n = nrow(aptos), y = log(aptos$delta))

fit_normal <- sampling(
  modelo_normal, data = datos_normal,
  chains = 4, iter = 2000, warmup = 1000,
  seed = 1234, refresh = 0
)

## ---- normal-resumen ----
resumen_normal <- resumen_stan(fit_normal, c("mu", "sigma"))
resumen_normal

## ---- guardar-normal-stan ----
guardar(resumen_normal, "resumen_normal")

## ---- normal-conjunta ----
draws_normal <- as_draws_df(fit_normal)
grafico_conjunta <- mcmc_scatter(fit_normal, pars = c("mu", "sigma"), alpha = 0.2) +
  mcmc_hist(fit_normal, pars = c("mu", "sigma")) +
  plot_layout(widths = c(1, 1.6))
guardar_figura(grafico_conjunta, "grafico_conjunta", alto = 5)
guardar_valor("cor_mu_sigma", cor(draws_normal$mu, draws_normal$sigma))

## ---- normal-funciones ----
tabla_funciones <- draws_normal %>%
  mutate(mediana_delta = exp(mu),
         media_delta = exp(mu + sigma^2 / 2)) %>%
  select(mediana_delta, media_delta, delta_nuevo) %>%
  summarise_draws(mean, sd, ~ quantile(.x, c(0.025, 0.975)))
tabla_funciones

## ---- guardar-normal-funciones ----
guardar(tabla_funciones, "tabla_funciones")

## ---- normal-ppc ----
set.seed(1)
idx <- sample(nrow(draws_normal), 50)
y_rep <- t(sapply(idx, function(s)
  rnorm(nrow(aptos), draws_normal$mu[s], draws_normal$sigma[s])))

grafico_ppc <- ppc_dens_overlay(log(aptos$delta), y_rep) +
  labs(x = "log(personas por vivienda)")
guardar_figura(grafico_ppc, "grafico_ppc", alto = 4.6)

# =============================================================================
# Multinomial-Dirichlet
# =============================================================================

## ---- multinomial-stan ----
conteos <- table(segmentos$clasificacion_cobertura)

modelo_multi <- stan_model("stan/multinomial_dirichlet.stan")
fit_multi <- sampling(
  modelo_multi,
  data = list(K = length(conteos), y = as.integer(conteos),
              alpha = rep(1, length(conteos))),
  chains = 4, iter = 2000, seed = 1234, refresh = 0
)

## ---- multinomial-resumen ----
alpha_post <- 1 + as.integer(conteos)
resumen_multi <- resumen_stan(fit_multi, "theta") %>%
  mutate(categoria = names(conteos),
         media_exacta = alpha_post / sum(alpha_post),
         desv_exacta = sqrt(media_exacta * (1 - media_exacta) / (sum(alpha_post) + 1)),
         .after = parametro) %>%
  select(parametro, categoria, media, media_exacta, desv_est, desv_exacta, Rhat)
resumen_multi

## ---- guardar-multinomial ----
guardar(resumen_multi, "resumen_multi")

## ---- multinomial-cor ----
theta <- as.data.frame(fit_multi) %>%
  select(starts_with("theta")) %>%
  setNames(names(conteos))
grafico_pares <- mcmc_pairs(theta, off_diag_args = list(size = 0.5, alpha = 0.2))
guardar_figura(grafico_pares, "grafico_pares", ancho = 10, alto = 6)
cor_theta <- cor(theta)
guardar(as.data.frame(cor_theta) %>% rownames_to_column("categoria"), "cor_theta")

# =============================================================================
# MCMC: Metropolis programado a mano
# =============================================================================

## ---- metropolis ----
log_posterior <- function(mu, log_sigma) {
  sigma <- exp(log_sigma)
  sum(dnorm(y, mu, sigma, log = TRUE)) +             # verosimilitud
    dnorm(mu, 0, 10, log = TRUE) +                    # previa de mu
    log(2) + dcauchy(sigma, 0, 2.5, log = TRUE) + log_sigma       # previa de sigma + jacobiano
}

metropolis <- function(inicio, n_iter = 3000, paso = 0.01) {
  cadena <- matrix(NA, n_iter, 2, dimnames = list(NULL, c("mu", "log_sigma")))
  accepted <- logical(n_iter)
  actual <- inicio
  lp_actual <- log_posterior(actual[1], actual[2])
  for (t in 1:n_iter) {
    propuesta <- actual + rnorm(2, 0, paso)
    lp_prop <- log_posterior(propuesta[1], propuesta[2])
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp_actual) {
      actual <- propuesta; lp_actual <- lp_prop; accepted[t] <- TRUE
    }
    cadena[t, ] <- actual
  }
  attr(cadena, "accepted") <- accepted
  cadena
}

## ---- metropolis-cadenas ----
set.seed(2026)
inicios <- list(c(-0.5, log(0.05)), c(2, log(1)), c(0, log(0.5)), c(1.5, log(0.1)))
cadenas <- map(inicios, metropolis)

grafico_cadenas <- map_dfr(seq_along(cadenas), ~ tibble(cadena = factor(.x), iteracion = 1:3000,
                                                        mu = cadenas[[.x]][, "mu"],
                                                        sigma = exp(cadenas[[.x]][, "log_sigma"]))) %>%
  pivot_longer(c(mu, sigma)) %>%
  ggplot(aes(iteracion, value, colour = cadena)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~ name, scales = "free_y") +
  labs(x = "Iteración", y = NULL, colour = "Cadena")
guardar_figura(grafico_cadenas, "grafico_cadenas", alto = 4.8)

## ---- metropolis-diagnostico ----
post <- map(cadenas, ~ .x[1001:3000, ])
arreglo <- simplify2array(post) %>% aperm(c(1, 3, 2))   # iteracion x cadena x parametro

diagnostico_metropolis <- tibble(
  parametro = c("mu", "log_sigma"),
  Rhat = c(rhat(arreglo[, , 1]), rhat(arreglo[, , 2])),
  ESS_bulk = round(c(ess_bulk(arreglo[, , 1]), ess_bulk(arreglo[, , 2]))),
  ESS_tail = round(c(posterior::ess_tail(arreglo[, , 1]), posterior::ess_tail(arreglo[, , 2]))),
  tasa_aceptacion = mean(unlist(map(cadenas, ~ attr(.x, "accepted")[1001:3000])))
)
diagnostico_metropolis

## ---- guardar-metropolis ----
guardar(diagnostico_metropolis, "diagnostico_metropolis")

## ---- autocorrelacion ----
acf_datos <- function(x, nombre) {
  a <- acf(x, lag.max = 60, plot = FALSE)
  tibble(algoritmo = nombre, rezago = drop(a$lag), autocorrelacion = drop(a$acf))
}
mu_stan <- as.array(fit_normal)[, 1, "mu"]
grafico_acf <- bind_rows(acf_datos(arreglo[, 1, 1], "Metropolis"),
                         acf_datos(mu_stan, "Stan (NUTS)")) %>%
  ggplot(aes(rezago, autocorrelacion, fill = algoritmo)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#c55a11", "#1f4e79")) +
  labs(x = "Rezago", y = "Autocorrelación de mu", fill = NULL) +
  theme(legend.position = "bottom")
guardar_figura(grafico_acf, "grafico_acf", alto = 4.5)

## ---- tamano-paso ----
set.seed(7)
experimento_paso <- map_dfr(c(0.001, 0.005, 0.01, 0.03, 0.1), function(p) {
  cs <- map(inicios, metropolis, paso = p)
  arr <- simplify2array(map(cs, ~ .x[1001:3000, ])) %>% aperm(c(1, 3, 2))
  tibble(paso = p, tasa_aceptacion = mean(unlist(map(cs, ~ attr(.x, "accepted")[1001:3000]))),
         Rhat_mu = rhat(arr[, , 1]), ESS_mu = round(ess_bulk(arr[, , 1])))
})
experimento_paso

## ---- guardar-tamano-paso ----
guardar(experimento_paso, "experimento_paso")

## ---- comparacion-ess ----
ess_stan <- summarise_draws(as_draws_array(fit_normal), "ess_bulk") %>%
  filter(variable %in% c("mu", "sigma"))

comparacion_ess <- tibble(
  algoritmo = c("Metropolis (4 x 2.000 draws)", "Stan NUTS (4 x 1.000 draws)"),
  ESS_mu = round(c(ess_bulk(arreglo[, , 1]), ess_stan$ess_bulk[match("mu", ess_stan$variable)])),
  ESS_sigma = round(c(ess_bulk(exp(arreglo[, , 2])), ess_stan$ess_bulk[match("sigma", ess_stan$variable)])))
comparacion_ess

## ---- guardar-comparacion-ess ----
guardar(comparacion_ess, "comparacion_ess")

## ---- diagnostic-stan ----
diagnostic_stan <- posterior::summarise_draws(
  posterior::as_draws_array(fit_normal), "rhat", "ess_bulk", "ess_tail", "mcse_mean"
) %>% filter(variable %in% c("mu", "sigma"))
guardar(diagnostic_stan, "diagnostico_stan")
nuts <- rstan::get_sampler_params(fit_normal, inc_warmup = FALSE)
diagnostic_nuts <- map_dfr(seq_along(nuts), function(i) {
  tibble(chain = i, divergences = sum(nuts[[i]][, "divergent__"]),
         max_treedepth_hits = sum(nuts[[i]][, "treedepth__"] >= 10))
})
guardar(diagnostic_nuts, "diagnostico_nuts")
writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")


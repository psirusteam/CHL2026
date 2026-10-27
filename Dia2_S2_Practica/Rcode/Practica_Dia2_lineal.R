##############################################################################
## Dia 2 - Sesion 2: Practica computacional en RStan                        ##
## Regresion bayesiana de la densidad: regiones, LOO, sigma por zona,       ##
## prediccion y chequeo predictivo previo.                                  ##
## Ejecutar con el directorio de trabajo en Dia2_S2_Practica/:              ##
##   source("Rcode/Practica_Dia2_lineal.R")                                 ##
##############################################################################

source("Rcode/utilidades.R")

## ---- paquetes ----
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_minimal(base_size = 15))
color_scheme_set("blue")

## ---- datos ----
magallanes <- readRDS("data/segmentos_magallanes.rds")
aysen <- readRDS("data/segmentos_aysen.rds")

covariables <- grep("_std$", names(magallanes), value = TRUE)

construir_X <- function(datos) {
  cbind(intercepto  = 1,
        zona_urbana = as.numeric(datos$zona == "urbana"),
        zona_aldea  = as.numeric(datos$zona == "aldea"),
        as.matrix(datos[, covariables]))
}

modelo_lineal <- stan_model("stan/regresion_lineal.stan")

# =============================================================================
# Ejercicio 1: efectos en dos regiones
# =============================================================================

## ---- ej1 ----
ajustar_lineal <- function(datos, seed) {
  aptos <- datos %>% filter(apto_modelo)
  X <- construir_X(aptos)
  stopifnot(all(is.finite(X)), all(is.finite(aptos$delta)), all(aptos$delta > 0),
            nrow(X) > ncol(X), qr(X)$rank == ncol(X))
  fit <- sampling(modelo_lineal,
                  data = list(n = nrow(X), p = ncol(X), X = X, y = log(aptos$delta), sd_beta = 1),
                  chains = 4, iter = 2000, warmup = 1000, seed = seed, refresh = 0, pars = c("beta", "sigma"))
  fit
}
fit_mag <- ajustar_lineal(magallanes, seed = 201)
fit_ays <- ajustar_lineal(aysen, seed = 202)
names_beta <- colnames(construir_X(filter(aysen, apto_modelo)))
region_draws <- function(fit) {
  beta <- as.data.frame(fit, pars = "beta")
  names(beta) <- names_beta
  beta$sigma <- as.matrix(fit, pars = "sigma")[, 1]
  beta$lp__ <- as.matrix(fit, pars = "lp__")[, 1]
  beta
}
draws_mag <- region_draws(fit_mag)
draws_ays <- region_draws(fit_ays)

## ---- ej1-grafico ----
efectos <- bind_rows(Magallanes = draws_mag, "Aysén" = draws_ays, .id = "region") %>%
  select(-intercepto, -sigma, -lp__) %>%
  pivot_longer(-region, names_to = "termino") %>%
  group_by(region, termino) %>%
  summarise(media = mean(value), li = quantile(value, 0.025), ls = quantile(value, 0.975), .groups = "drop")

grafico_ej1 <- ggplot(efectos, aes(media, str_remove(termino, "_std$"), colour = region)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = li, xmax = ls), position = position_dodge(width = 0.6)) +
  labs(x = "Coeficiente", y = NULL, colour = "Región")

## ---- guardar-ej1 ----
guardar_figura(grafico_ej1, "grafico_ej1", alto = 6)
signos <- efectos %>%
  select(region, termino, media) %>%
  pivot_wider(names_from = region, values_from = media) %>%
  mutate(termino = str_remove(termino, "_std$"),
         mismo_signo = sign(`Aysén`) == sign(Magallanes))
guardar(signos, "signos")
guardar(tibble(region = c("Aysén", "Magallanes"),
               sigma = c(mean(draws_ays$sigma), mean(draws_mag$sigma))), "sigmas_region")

# =============================================================================
# Ejercicio 2: comparacion de modelos con LOO
# =============================================================================

## ---- ej2 ----
aptos_ays <- aysen %>% filter(apto_modelo)
y <- log(aptos_ays$delta)
X_B <- construir_X(aptos_ays)
X_A <- X_B[, c("intercepto", "zona_urbana", "zona_aldea")]

ajustar_X <- function(X, seed) {
  sampling(modelo_lineal, data = list(n = nrow(X), p = ncol(X), X = X, y = y, sd_beta = 1),
           chains = 4, iter = 2000, warmup = 1000, seed = seed, refresh = 0, pars = c("beta", "sigma"))
}
fit_A <- ajustar_X(X_A, seed = 203)
fit_B <- fit_ays  # same data, design and prior: reuse the regional fit

## ---- ej2-loo ----
log_lik_lineal <- function(fit, X) {
  beta <- as.matrix(fit, pars = "beta"); sigma <- as.matrix(fit, pars = "sigma")[, 1]
  mu <- beta %*% t(X)                                       # draws x segmentos
  dnorm(matrix(y, nrow(mu), ncol(mu), byrow = TRUE), mu, sigma, log = TRUE)
}
compute_loo <- function(log_lik, fit) {
  dims <- dim(as.array(fit, pars = "beta"))
  chain_id <- rep(seq_len(dims[2]), each = dims[1])
  scaled_likelihood <- exp(sweep(log_lik, 2, apply(log_lik, 2, max), "-"))
  r_eff <- loo::relative_eff(scaled_likelihood, chain_id = chain_id, cores = 2)
  loo::loo(log_lik, r_eff = r_eff, cores = 2)
}
loo_A <- compute_loo(log_lik_lineal(fit_A, X_A), fit_A)
loo_B <- compute_loo(log_lik_lineal(fit_B, X_B), fit_B)
comparacion_loo <- loo_compare(list(`A: zona` = loo_A, `B: zona + satelitales` = loo_B))
comparacion_loo

## ---- guardar-ej2 ----
guardar_texto(comparacion_loo, "comparacion_loo")
guardar(tibble(modelo = c("A: zona", "B: zona + satelitales"),
               sigma = c(mean(as.matrix(fit_A, pars = "sigma")), mean(as.matrix(fit_B, pars = "sigma"))),
               elpd_loo = c(loo_A$estimates["elpd_loo", 1], loo_B$estimates["elpd_loo", 1]),
               p_loo = c(loo_A$estimates["p_loo", 1], loo_B$estimates["p_loo", 1]),
               k_mayor_0.7 = c(sum(loo_A$diagnostics$pareto_k > 0.7), sum(loo_B$diagnostics$pareto_k > 0.7))),
        "tabla_loo")
delta_pointwise <- loo_B$pointwise[, "elpd_loo"] - loo_A$pointwise[, "elpd_loo"]
guardar_valor("elpd_diff", sum(delta_pointwise))
guardar_valor("se_diff", sqrt(length(delta_pointwise) * var(delta_pointwise)))

# =============================================================================
# Ejercicio 3: una sigma por zona
# =============================================================================

## ---- ej3 ----
zonas <- c("rural", "aldea", "urbana")
modelo_sz <- stan_model("stan/regresion_sigma_zona.stan")
fit_sz <- sampling(
  modelo_sz,
  data = list(n = nrow(X_B), p = ncol(X_B), X = X_B, y = y,
              J = 3, zona = match(aptos_ays$zona, zonas)),
  chains = 4, iter = 2000, warmup = 1000, seed = 204, refresh = 0
)

## ---- ej3-resultado ----
sigmas_zona <- summarise_draws(as_draws_df(fit_sz), "mean", ~ quantile(.x, c(0.025, 0.975))) %>%
  filter(str_starts(variable, "sigma")) %>%
  mutate(zona = zonas, .after = variable)
loo_sz <- compute_loo(loo::extract_log_lik(fit_sz), fit_sz)
comparacion_sz <- loo_compare(list(`B: sigma común` = loo_B, `C: sigma por zona` = loo_sz))
sigmas_zona
comparacion_sz

## ---- guardar-ej3 ----
guardar(sigmas_zona, "sigmas_zona")
guardar_texto(comparacion_sz, "comparacion_sz")

## ---- ej3-ppc ----
y_rep_sz <- as.matrix(fit_sz, pars = "y_rep")
grafico_ej3 <- ppc_stat_grouped(y, y_rep_sz, group = aptos_ays$zona, stat = "sd") +
  labs(subtitle = "Modelo con una sigma por zona")
guardar_figura(grafico_ej3, "grafico_ej3", alto = 4.6)

# =============================================================================
# Ejercicio 4: densidad y poblacion en segmentos de cobertura incompleta
# =============================================================================

## ---- ej4 ----
beta_B <- as.matrix(fit_B, pars = "beta"); sigma_B <- as.matrix(fit_B, pars = "sigma")[, 1]
X_todos <- construir_X(aysen)
stopifnot(identical(colnames(X_todos), colnames(X_B)), all(is.finite(X_todos)),
          !anyNA(aysen$n_vpo), all(aysen$n_vpo >= 0))

set.seed(4)
eta <- beta_B %*% t(X_todos)
delta_pred <- exp(eta + matrix(rnorm(length(eta)), nrow(eta)) * sigma_B)
y_pred <- delta_pred * matrix(aysen$n_vpo, nrow(eta), ncol(eta), byrow = TRUE)
stopifnot(all(is.finite(y_pred)))

ej4 <- map_dfr(levels(aysen$clasificacion_cobertura), function(k) {
  cols <- aysen$clasificacion_cobertura == k
  total <- rowSums(y_pred[, cols, drop = FALSE])
  tibble(cobertura = k, censado = sum(aysen$n_per[cols]),
         predicho = mean(total), q025 = quantile(total, 0.025), q975 = quantile(total, 0.975),
         V_cero = sum(aysen$n_vpo[cols] == 0))
})
ej4

## ---- guardar-ej4 ----
guardar(ej4, "ej4")

# =============================================================================
# Ejercicio 5: chequeo predictivo previo
# =============================================================================

## ---- ej5 ----
simular_previa <- function(sd_beta, media_b0, sd_b0, sim_sigma, seed, n_sim = 200) {
  set.seed(seed)
  map_dfr(1:n_sim, function(s) {
    beta <- c(rnorm(1, media_b0, sd_b0), rnorm(ncol(X_todos) - 1, 0, sd_beta))
    sigma <- sim_sigma()
    tibble(log_delta = drop(X_todos %*% beta) + rnorm(nrow(X_todos), 0, sigma))
  })
}
previa_vaga <- simular_previa(1, 0, 1, function() abs(rcauchy(1, 0, 2.5)), seed = 501)
previa_info <- simular_previa(0.2, log(2.5), 0.5, function() abs(rnorm(1, 0, 0.5)), seed = 502)

ej5 <- bind_rows(`β ~ N(0, 1), σ ~ Cauchy⁺` = previa_vaga, `Previas informativas` = previa_info, .id = "previa") %>%
  group_by(previa) %>%
  summarise(`% δ > 20` = 100 * mean(log_delta > log(20)), `% δ < 1` = 100 * mean(log_delta < 0),
            `mediana δ` = exp(median(log_delta)), `P90 δ` = exp(quantile(log_delta, 0.90)),
            pct_overflow = 100 * mean(log_delta > log(.Machine$double.xmax)),
            pct_plot_omitted = 100 * mean(log_delta <= log(0.01) | log_delta >= log(100)))
ej5

## ---- guardar-ej5 ----
guardar(ej5, "ej5")
grafico_ej5 <- bind_rows(`β ~ N(0, 1), σ ~ Cauchy⁺` = previa_vaga, `Previas informativas` = previa_info, .id = "previa") %>%
  filter(log_delta > log(0.01), log_delta < log(100)) %>%
  mutate(delta = exp(log_delta)) %>%
  ggplot(aes(delta, fill = previa)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  scale_x_log10() +
  geom_vline(xintercept = c(1, 6), linetype = "dashed") +
  labs(x = "δ simulada desde la previa (escala log)", y = "Frecuencia", fill = NULL,
       subtitle = "Líneas: referencias ilustrativas 1 y 6; se muestran solo valores entre 0,01 y 100") +
  theme(legend.position = "bottom")
guardar_figura(grafico_ej5, "grafico_ej5", alto = 4.5)

## ---- regional-differences ----
regional_differences <- purrr::map_dfr(setdiff(names_beta, "intercepto"), function(term) {
  difference <- draws_mag[[term]] - draws_ays[[term]]
  tibble(term = term, mean_difference = mean(difference),
         q025 = quantile(difference, 0.025), q975 = quantile(difference, 0.975),
         probability_positive = mean(difference > 0))
})
guardar(regional_differences, "diferencias_region")

## ---- diagnostics ----
fits <- list(magallanes = fit_mag, aysen_B = fit_B, aysen_A = fit_A, aysen_C = fit_sz)
diagnostics <- purrr::imap_dfr(fits, function(fit, name) {
  posterior::summarise_draws(posterior::as_draws_array(as.array(fit, pars = c("beta", "sigma"))),
                           "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
    mutate(fit = name, .before = 1)
})
guardar(diagnostics, "diagnostico_ajustes")
nuts <- purrr::imap_dfr(fits, function(fit, name) {
  chains <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  tibble(fit = name, divergences = sum(vapply(chains, function(x) sum(x[, "divergent__"]), numeric(1))),
         max_treedepth_hits = sum(vapply(chains, function(x) sum(x[, "treedepth__"] >= 10), numeric(1))))
})
guardar(nuts, "diagnostico_nuts")
guardar(list(A = loo_A, B = loo_B, C = loo_sz), "loo_models")
loo_diagnostics <- purrr::imap_dfr(list(A = loo_A, B = loo_B, C = loo_sz), function(x, name) {
  threshold <- min(1 - 1/log10(attr(x, "dims")[1]), 0.7)
  tibble(model = name, k_threshold = threshold, k_max = max(loo::pareto_k_values(x)),
         n_flagged = length(loo::pareto_k_ids(x)),
         n_k_above_one = sum(loo::pareto_k_values(x) > 1))
})
guardar(loo_diagnostics, "diagnostico_loo")
writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")



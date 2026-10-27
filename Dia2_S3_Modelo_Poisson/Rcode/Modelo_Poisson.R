##############################################################################
## Dia 2 - Sesion 3: Modelo Poisson con vinculo loglineal                  ##
## Poisson con offset, sobredispersion y Poisson-lognormal en Aysen.       ##
## Ejecutar con el directorio de trabajo en Dia2_S3_Modelo_Poisson/:       ##
##   source("Rcode/Modelo_Poisson.R")                                       ##
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
stopifnot(all(is.finite(aptos$n_vpo)), all(aptos$n_vpo > 0),
          all(is.finite(aptos$n_per)), all(aptos$n_per > 0),
          all(aptos$n_per == floor(aptos$n_per)))

## ---- offset ----
ajuste_offset <- glm(n_per ~ log(n_vpo), family = poisson, data = aptos)
ajuste_offset_quasi <- glm(n_per ~ log(n_vpo), family = quasipoisson, data = aptos)
tabla_offset <- bind_rows(Poisson = broom::tidy(ajuste_offset),
                          quasiPoisson = broom::tidy(ajuste_offset_quasi), .id = "model") %>%
  select(model, term, estimate, std.error)
tabla_offset

## ---- guardar-offset ----
guardar(tabla_offset, "tabla_offset")
grafico_offset <- ggplot(aptos, aes(n_vpo, n_per)) +
  geom_point(alpha = 0.3, colour = "#2e75b6") +
  geom_abline(intercept = log10(sum(aptos$n_per) / sum(aptos$n_vpo)), slope = 1, colour = "#c00000", linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Viviendas ocupadas (escala log)", y = "Personas (escala log)",
       subtitle = "Línea roja: pendiente 1 (personas proporcionales a viviendas)")
guardar_figura(grafico_offset, "grafico_offset", alto = 4.2)

## ---- x-glm ----
covariables <- grep("_std$", names(segmentos), value = TRUE)
construir_X <- function(datos) {
  cbind(intercepto  = 1,
        zona_urbana = as.numeric(datos$zona == "urbana"),
        zona_aldea  = as.numeric(datos$zona == "aldea"),
        as.matrix(datos[, covariables]))
}
X <- construir_X(aptos)
stopifnot(all(is.finite(X)), nrow(X) > ncol(X), qr(X)$rank == ncol(X))

ajuste_glm <- glm(aptos$n_per ~ X - 1 + offset(log(aptos$n_vpo)), family = poisson)
dispersion_glm <- sum(residuals(ajuste_glm, type = "pearson")^2) / ajuste_glm$df.residual
dispersion_glm

## ---- guardar-glm ----
guardar_valor("dispersion_glm", dispersion_glm)

## ---- ajuste-poisson ----
modelo_poisson <- stan_model("stan/poisson_offset.stan")

datos_stan <- list(n = nrow(X), p = ncol(X), X = X,
                   y = aptos$n_per, log_v = log(aptos$n_vpo))

set.seed(303)
init_poisson <- lapply(seq_len(4), function(chain) {
  list(beta = unname(coef(ajuste_glm) + rnorm(ncol(X), 0, 0.5 * sqrt(diag(vcov(ajuste_glm))))))
})
fit_poisson <- sampling(modelo_poisson, data = datos_stan, init = init_poisson,
                        chains = 4, iter = 2000, warmup = 1000, seed = 301, refresh = 0,
                        control = list(metric = "dense_e", adapt_delta = 0.95, max_treedepth = 12))

saveRDS(fit_poisson, "resultados/fit_poisson.rds")

## ---- resultado-poisson ----
resultado_poisson <- summarise_draws(as_draws_array(as.array(fit_poisson, pars = "beta")),
                  "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
  filter(str_starts(variable, "beta")) %>%
  mutate(termino = str_remove(colnames(X), "_std$"), glm = coef(ajuste_glm),
         error_glm = sqrt(diag(vcov(ajuste_glm))), .after = variable)
resultado_poisson

## ---- guardar-poisson ----
guardar(resultado_poisson, "resultado_poisson")

## ---- ppc-poisson ----
dispersion <- function(y, mu) sum((y - mu)^2 / mu) / (length(y) - ncol(X))

mu_hat <- colMeans(exp(as.matrix(fit_poisson, pars = "beta") %*% t(X))) * aptos$n_vpo
y_rep <- as.matrix(fit_poisson, pars = "y_rep")
stopifnot(all(as.matrix(fit_poisson, pars = "rng_out_of_range") == 0), all(y_rep >= 0))
set.seed(1)
phi_rep <- apply(y_rep[sample(nrow(y_rep), 500), ], 1, dispersion, mu = mu_hat)
phi_obs <- dispersion(aptos$n_per, mu_hat)

## ---- guardar-ppc-poisson ----
guardar_valor("phi_obs", phi_obs)
guardar_valor("phi_rep_q95", unname(quantile(phi_rep, 0.95)))
grafico_ppc_poisson <- ggplot(tibble(phi_rep), aes(phi_rep)) +
  geom_histogram(bins = 40, fill = "#9ecae1", colour = "white") +
  geom_vline(xintercept = phi_obs, colour = "#c00000", linewidth = 1.2) +
  labs(x = "Índice de dispersión de Pearson", y = "Réplicas",
       subtitle = "Histograma: réplicas del modelo Poisson · Línea roja: datos observados")
guardar_figura(grafico_ppc_poisson, "grafico_ppc_poisson", alto = 4.2)

## ---- ajuste-pln ----
modelo_pln <- stan_model("stan/poisson_lognormal.stan")

set.seed(304)
initial_lm <- lm(log(aptos$n_per / aptos$n_vpo) ~ X - 1)
initial_beta <- unname(coef(initial_lm))
initial_scales <- c(0.08, 0.10, 0.12, 0.15)
init_pln <- lapply(seq_len(4), function(chain) {
  beta <- initial_beta + rnorm(ncol(X), 0, 0.1 * sqrt(diag(vcov(initial_lm))))
  sigma <- initial_scales[chain]
  z <- drop((log(aptos$n_per / aptos$n_vpo) - X %*% beta) / sigma)
  list(beta = beta, sigma = sigma, z = z + rnorm(nrow(X), 0, 0.05))
})
fit_pln <- sampling(modelo_pln, data = datos_stan, init = init_pln,
                    chains = 4, iter = 2000, warmup = 1000, seed = 302, refresh = 0,
                    control = list(adapt_delta = 0.95, max_treedepth = 12),
                    pars = c("beta", "sigma", "z", "log_delta", "y_rep", "rng_out_of_range"))

saveRDS(fit_pln, "resultados/fit_pln.rds")

## ---- resultado-pln ----
sampler_pln <- rstan::get_sampler_params(fit_pln, inc_warmup = FALSE)
diagnostico_pln <- tibble(
  divergencias = sum(vapply(sampler_pln, function(x) sum(x[, "divergent__"]), numeric(1))),
  max_treedepth = sum(vapply(sampler_pln, function(x) sum(x[, "treedepth__"] >= 12), numeric(1))))
resumen_pln <- summarise_draws(as_draws_array(as.array(fit_pln, pars = c("beta", "sigma"))),
                  "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
  filter(variable == "sigma" | str_starts(variable, "beta\\[[1-3]\\]"))
diagnostico_pln
resumen_pln

## ---- guardar-pln ----
guardar(diagnostico_pln, "diagnostico_pln")
guardar(resumen_pln, "resumen_pln")
sigma_pln <- mean(as.matrix(fit_pln, pars = "sigma"))
guardar_valor("sigma_pln", sigma_pln)
guardar_valor("sigma_lineal", sigma(lm(log(aptos$delta) ~ X - 1)))

## ---- momentos-pln ----
# Momentos marginales del Poisson-lognormal: E(Y) = V exp(eta + sigma^2/2),
# Var(Y) = E(Y) + E(Y)^2 (exp(sigma^2) - 1). Se comparan con los de las replicas.
beta_pln <- as.matrix(fit_pln, pars = "beta")
sigma_draws <- as.matrix(fit_pln, pars = "sigma")[, 1]
eta_draws <- beta_pln %*% t(X)
eta_hat <- matrix(colMeans(eta_draws), nrow = 1)
mean_draws <- exp(eta_draws + sigma_draws^2 / 2) *
  matrix(aptos$n_vpo, nrow(eta_draws), ncol(eta_draws), byrow = TRUE)
variance_draws <- mean_draws + mean_draws^2 * expm1(sigma_draws^2)
media_teorica <- colMeans(mean_draws)
var_teorica <- colMeans(variance_draws) + apply(mean_draws, 2, var)
y_rep_pln <- as.matrix(fit_pln, pars = "y_rep")
stopifnot(all(as.matrix(fit_pln, pars = "rng_out_of_range") == 0),
          all(y_rep_pln >= 0), all(is.finite(var_teorica)))
tabla_momentos <- tibble(
  cantidad = c("Media de Y (promedio de segmentos)", "Varianza de Y (promedio)", "Varianza / media"),
  teorica = c(mean(media_teorica), mean(var_teorica), mean(var_teorica / media_teorica)),
  replicas = c(mean(colMeans(y_rep_pln)), mean(apply(y_rep_pln, 2, var)),
               mean(apply(y_rep_pln, 2, var) / colMeans(y_rep_pln))))
tabla_momentos

## ---- guardar-momentos ----
guardar(tabla_momentos, "tabla_momentos")

## ---- ppc-pln ----
set.seed(2)
phi_rep_pln <- apply(y_rep_pln[sample(nrow(y_rep_pln), 500), ], 1, dispersion, mu = mu_hat)

grafico_ppc_pln <- bind_rows(tibble(modelo = "Poisson", phi = phi_rep),
                             tibble(modelo = "Poisson–lognormal", phi = phi_rep_pln)) %>%
  ggplot(aes(phi, fill = modelo)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity", colour = "white") +
  geom_vline(xintercept = phi_obs, colour = "#c00000", linewidth = 1.2) +
  scale_fill_manual(values = c("#9ecae1", "#548235")) +
  labs(x = "Índice de dispersión de Pearson", y = "Réplicas", fill = NULL,
       subtitle = "Línea roja: valor observado") +
  theme(legend.position = "bottom")
guardar_figura(grafico_ppc_pln, "grafico_ppc_pln", alto = 4.5)

## ---- ppc-densidad ----
set.seed(1)
idx <- sample(nrow(y_rep_pln), 50)
v <- matrix(aptos$n_vpo, length(idx), nrow(aptos), byrow = TRUE)

grafico_ppc_densidad <- (ppc_dens_overlay(aptos$n_per / aptos$n_vpo, y_rep_pln[idx, ] / v) +
                           labs(subtitle = "Poisson–lognormal")) +
  (ppc_dens_overlay(aptos$n_per / aptos$n_vpo, y_rep[idx, ] / v) + labs(subtitle = "Poisson")) &
  coord_cartesian(xlim = c(1, 4.5)) & labs(x = "Personas por vivienda")
guardar_figura(grafico_ppc_densidad, "grafico_ppc_densidad", alto = 4.6)

## ---- encogimiento ----
# Densidad posterior de cada segmento frente a la observada y a la de las covariables.
# Aproximacion: log(Y/V) tiene varianza ~ 1/Y, asi que el peso del dato es
#   w_i = sigma^2 / (sigma^2 + 1/Y_i)
log_delta_draws <- as.matrix(fit_pln, pars = "log_delta")
log_delta_post <- colMeans(log_delta_draws)
delta_median <- exp(apply(log_delta_draws, 2, median))
eta_cov <- eta_hat[1, ]
w_aprox <- sigma_pln^2 / (sigma_pln^2 + 1 / aptos$n_per)
w_stan <- (log_delta_post - eta_cov) / (log(aptos$n_per / aptos$n_vpo) - eta_cov)

tabla_encogimiento <- tibble(viviendas = aptos$n_vpo, w_aprox, w_stan) %>%
  filter(abs(log(aptos$n_per / aptos$n_vpo) - eta_cov) > 0.05) %>%     # evita divisiones inestables
  mutate(tamano = cut(viviendas, c(0, 10, 20, 40, Inf), labels = c("1–10", "11–20", "21–40", "más de 40"))) %>%
  group_by(tamano) %>%
  summarise(segmentos = n(), peso_dato_aprox = median(w_aprox), ratio_implicito_stan = median(w_stan))
tabla_encogimiento

## ---- guardar-encogimiento ----
guardar(tabla_encogimiento, "tabla_encogimiento")
grafico_encogimiento <- tibble(observada = aptos$n_per / aptos$n_vpo, posterior = delta_median,
                               viviendas = aptos$n_vpo) %>%
  ggplot(aes(observada, posterior, colour = log10(viviendas))) +
  geom_abline(linetype = "dashed") +
  geom_point(alpha = 0.6) +
  scale_colour_viridis_c(labels = function(x) round(10^x), name = "Viviendas") +
  labs(x = "Densidad observada Y/V", y = "Densidad posterior (mediana)")
guardar_figura(grafico_encogimiento, "grafico_encogimiento", alto = 4.8)

## ---- diagnostics ----
diagnostic_draws <- list(
  Poisson = as.array(fit_poisson, pars = "beta"),
  Poisson_lognormal = as.array(fit_pln, pars = c("beta", "sigma", "z")))
diagnostics <- purrr::imap_dfr(diagnostic_draws, function(draws, name) {
  posterior::summarise_draws(posterior::as_draws_array(draws),
                           "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
    mutate(model = name, .before = 1)
})
guardar(diagnostics, "diagnostico_completo")
nuts <- purrr::imap_dfr(list(Poisson = fit_poisson, Poisson_lognormal = fit_pln), function(fit, name) {
  chains <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  tibble(model = name,
         divergences = sum(vapply(chains, function(x) sum(x[, "divergent__"]), numeric(1))),
         max_treedepth_hits = sum(vapply(chains, function(x) sum(x[, "treedepth__"] >= 12), numeric(1))),
         invalid_replications = sum(as.matrix(fit, pars = "rng_out_of_range")))
})
guardar(nuts, "diagnostico_nuts")
writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")


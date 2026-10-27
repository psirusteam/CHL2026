##############################################################################
## Dia 1 - Sesion 4: Practica computacional en RStan                        ##
## Fundamentos de la inferencia bayesiana con los segmentos de Magallanes.  ##
## Ejecutar con el directorio de trabajo en Dia1_S4_Practica/:              ##
##   source("Rcode/Practica_Dia1.R")                                        ##
##############################################################################

source("Rcode/utilidades.R")

## ---- paquetes ----
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_minimal(base_size = 15))
color_scheme_set("blue")

## ---- datos ----
magallanes <- readRDS("data/segmentos_magallanes.rds")
aysen <- readRDS("data/segmentos_aysen.rds")

glimpse(select(magallanes, edid, zona, n_per, n_vpo, n_vc, clasificacion_cobertura,
               delta, apto_modelo))

## ---- guardar-datos ----
guardar_texto(glimpse(select(magallanes, edid, zona, n_per, n_vpo, n_vc, clasificacion_cobertura,
                             delta, apto_modelo)), "glimpse_magallanes")

# =============================================================================
# Ejercicio 1: proporcion de segmentos con viviendas colectivas
# =============================================================================

## ---- ej1 ----
stopifnot(!anyNA(magallanes$n_vc), !anyNA(magallanes$zona))
y <- as.integer(magallanes$n_vc > 0)

posterior_beta <- function(a, b) {
  c(media = (a + sum(y)) / (a + b + length(y)),
    q025 = qbeta(0.025, a + sum(y), b + length(y) - sum(y)),
    q975 = qbeta(0.975, a + sum(y), b + length(y) - sum(y)))
}

modelo_bern <- stan_model("stan/bernoulli_beta.stan")
fit_ej1 <- sampling(modelo_bern, data = list(n = length(y), y = y, a = 1, b = 1),
                    chains = 4, iter = 2000, warmup = 1000, seed = 101, refresh = 0)

fit_ej1_sens <- sampling(modelo_bern,
  data = list(n = length(y), y = y, a = 1, b = 99),
  chains = 4, iter = 2000, warmup = 1000, seed = 102, refresh = 0)

## ---- ej1-resultado ----
summary_theta <- function(fit) {
  draws <- as.data.frame(fit)$theta
  c(media = mean(draws), q025 = unname(quantile(draws, 0.025)),
    q975 = unname(quantile(draws, 0.975)))
}
ej1 <- rbind(`Analítica Beta(1, 1)` = posterior_beta(1, 1),
             `Stan Beta(1, 1)` = summary_theta(fit_ej1),
             `Analítica Beta(1, 99)` = posterior_beta(1, 99),
             `Stan Beta(1, 99)` = summary_theta(fit_ej1_sens)) %>%
  as.data.frame() %>% rownames_to_column("metodo")
ej1

## ---- guardar-ej1 ----
guardar(ej1, "ej1")
guardar_valor("n_mag", length(y))
guardar_valor("colectivos_mag", sum(y))

# =============================================================================
# Ejercicio 2: personas por vivienda en dos regiones
# =============================================================================

## ---- ej2 ----
modelo_pg <- stan_model("stan/poisson_gamma.stan")

ajustar_pg <- function(datos, seed) {
  completos <- datos %>% filter(clasificacion_cobertura == "Completa", !excluido_colectivo, n_vpo > 0)
  sampling(modelo_pg,
           data = list(n = nrow(completos), y = completos$n_per, v = completos$n_vpo,
                       a = 0.01, b = 0.01),
           chains = 4, iter = 2000, warmup = 1000, seed = seed, refresh = 0)
}
fit_mag <- ajustar_pg(magallanes, seed = 201)
fit_ays <- ajustar_pg(aysen, seed = 202)

## ---- ej2-resultado ----
d_mag <- as.data.frame(fit_mag)$delta
d_ays <- as.data.frame(fit_ays)$delta
diferencia <- d_mag - d_ays
ej2 <- tibble(delta_magallanes = mean(d_mag), delta_aysen = mean(d_ays),
              dif_media = mean(diferencia),
              q025 = quantile(diferencia, 0.025), q975 = quantile(diferencia, 0.975),
              `Pr(dif > 0)` = mean(diferencia > 0))
ej2

## ---- guardar-ej2 ----
guardar(ej2, "ej2")
grafico_ej2 <- ggplot(tibble(diferencia), aes(diferencia)) +
  geom_histogram(bins = 50, fill = "#2e75b6", colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "δ Magallanes − δ Aysén (personas por vivienda)", y = "Draws")
guardar_figura(grafico_ej2, "grafico_ej2", alto = 4.2)

## ---- ej2-dispersion ----
completos_mag <- magallanes %>%
  filter(clasificacion_cobertura == "Completa", !excluido_colectivo, n_vpo > 0)
mu <- mean(d_mag) * completos_mag$n_vpo
dispersion_mag <- sum((completos_mag$n_per - mu)^2 / mu) / (nrow(completos_mag) - 1)
dispersion_mag

## ---- guardar-ej2-dispersion ----
guardar_valor("dispersion_mag", dispersion_mag)

# =============================================================================
# Ejercicio 3: densidad por zona
# =============================================================================

## ---- ej3 ----
modelo_normal <- stan_model("stan/normal_media_varianza.stan")

ajustar_zona <- function(z, seed) {
  datos_z <- magallanes %>% filter(apto_modelo, zona == z)
  stopifnot(nrow(datos_z) > 2, all(is.finite(datos_z$delta)), all(datos_z$delta > 0),
            var(log(datos_z$delta)) > 0)
  sampling(modelo_normal, data = list(n = nrow(datos_z), y = log(datos_z$delta)),
           chains = 4, iter = 2000, warmup = 1000, seed = seed, refresh = 0)
}
fit_urbana <- ajustar_zona("urbana", seed = 301)
fit_rural <- ajustar_zona("rural", seed = 302)

## ---- ej3-resultado ----
urb <- as.data.frame(fit_urbana); rur <- as.data.frame(fit_rural)
razon <- exp(urb$mu - rur$mu)

ej3 <- tibble(parametro = c("razón de medianas urbana/rural", "sigma urbana", "sigma rural",
                            "Pr(sigma rural > sigma urbana)"),
              media = c(mean(razon), mean(urb$sigma), mean(rur$sigma), mean(rur$sigma > urb$sigma)),
              q025 = c(quantile(razon, 0.025), quantile(urb$sigma, 0.025), quantile(rur$sigma, 0.025), NA),
              q975 = c(quantile(razon, 0.975), quantile(urb$sigma, 0.975), quantile(rur$sigma, 0.975), NA))
ej3

## ---- guardar-ej3 ----
guardar(ej3, "ej3")

## ---- ej3-ppc ----
y_urbana <- log(filter(magallanes, apto_modelo, zona == "urbana")$delta)
set.seed(1)
idx <- sample(nrow(urb), 50)
y_rep <- t(sapply(idx, function(s) rnorm(length(y_urbana), urb$mu[s], urb$sigma[s])))

grafico_ej3 <- ppc_dens_overlay(y_urbana, y_rep) +
  labs(x = "log(personas por vivienda), zona urbana")

## ---- guardar-ej3-ppc ----
guardar_figura(grafico_ej3, "grafico_ej3", alto = 4.6)

# =============================================================================
# Ejercicio 4: segmentos por zona
# =============================================================================

## ---- ej4 ----
zone_levels <- c("urbana", "aldea", "rural")
stopifnot(all(magallanes$zona %in% zone_levels))
conteos <- table(factor(magallanes$zona, levels = zone_levels))
modelo_multi <- stan_model("stan/multinomial_dirichlet.stan")
fit_ej4 <- sampling(modelo_multi,
                    data = list(K = length(conteos), y = as.integer(conteos), alpha = rep(1, length(conteos))),
                    chains = 4, iter = 2000, warmup = 1000, seed = 401, refresh = 0)

## ---- ej4-resultado ----
theta <- as.data.frame(fit_ej4) %>% select(starts_with("theta")) %>% setNames(names(conteos))

ej4 <- tibble(zona = names(conteos), segmentos = as.integer(conteos),
              media_stan = colMeans(theta),
              media_analitica = (1 + as.integer(conteos)) / (3 + sum(conteos)))
ej4
prob_urbana <- mean(theta$urbana > 0.75)
alpha_post <- 1 + as.integer(conteos)
urban_index <- match("urbana", names(conteos))
prob_urbana_exacta <- pbeta(0.75, alpha_post[urban_index],
                            sum(alpha_post) - alpha_post[urban_index], lower.tail = FALSE)

## ---- guardar-ej4 ----
guardar(ej4, "ej4")
guardar_valor("prob_urbana", prob_urbana)
guardar_valor("prob_urbana_exacta", prob_urbana_exacta)

# =============================================================================
# Ejercicio 5 (desafio): Metropolis para la densidad
# =============================================================================

## ---- ej5 ----
log_post <- function(log_d, y, v, a = 0.01, b = 0.01) {
  d <- exp(log_d)
  sum(dpois(y, d * v, log = TRUE)) + dgamma(d, a, b, log = TRUE) + log_d   # + jacobiano
}
metropolis_pg <- function(inicio, y, v, n_iter = 5000, paso = 0.01) {
  cadena <- numeric(n_iter); accepted <- logical(n_iter)
  actual <- inicio; lp <- log_post(actual, y, v)
  for (t in 1:n_iter) {
    prop <- actual + rnorm(1, 0, paso); lp_prop <- log_post(prop, y, v)
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) {
      actual <- prop; lp <- lp_prop; accepted[t] <- TRUE
    }
    cadena[t] <- actual
  }
  draws <- exp(cadena)
  attr(draws, "accepted") <- accepted
  draws
}

## ---- ej5-resultado ----
set.seed(5)
ej5 <- map_dfr(c(0.001, 0.01, 0.05), function(paso) {
  chains_full <- lapply(log(c(1, 2, 3, 4)), metropolis_pg, paso = paso,
                    y = completos_mag$n_per, v = completos_mag$n_vpo)
  cadenas <- do.call(cbind, lapply(chains_full, function(x) x[1001:5000]))
  tibble(paso = paso, media_mcmc = mean(cadenas),
         exacta = (0.01 + sum(completos_mag$n_per)) / (0.01 + sum(completos_mag$n_vpo)),
         tasa_aceptacion = mean(unlist(lapply(chains_full, function(x) attr(x, "accepted")[1001:5000]))),
         Rhat = posterior::rhat(cadenas), ESS_bulk = round(posterior::ess_bulk(cadenas)),
         ESS_tail = round(posterior::ess_tail(cadenas)))
})
ej5

## ---- guardar-ej5 ----
guardar(ej5, "ej5")

## ---- diagnostics ----
fits <- list(ej1_uniform = fit_ej1, ej1_sensitivity = fit_ej1_sens,
             ej2_magallanes = fit_mag, ej2_aysen = fit_ays,
             ej3_urban = fit_urbana, ej3_rural = fit_rural, ej4_zones = fit_ej4)
diagnostic_tables <- purrr::imap(fits, function(fit, name) {
  posterior::summarise_draws(posterior::as_draws_array(fit),
                           "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
    filter(grepl("^(theta|delta$|mu$|sigma$)", variable)) %>%
    mutate(fit = name, .before = 1)
})
guardar(bind_rows(diagnostic_tables), "diagnostico_ajustes")
nuts_diagnostics <- purrr::imap_dfr(fits, function(fit, name) {
  chains <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  tibble(fit = name, divergences = sum(vapply(chains, function(x) sum(x[, "divergent__"]), numeric(1))),
         max_treedepth_hits = sum(vapply(chains, function(x) sum(x[, "treedepth__"] >= 10), numeric(1))))
})
guardar(nuts_diagnostics, "diagnostico_nuts")
writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")

# D4S1: diagnostico del modelo de Magallanes, sin repetir Stan en clase.
# Ejecutar desde la carpeta de la sesion; las cadenas se preparan antes.

## ---- diagnostico-paquetes ----
library(posterior)
library(bayesplot)
library(ggplot2)
theme_set(theme_minimal(base_size = 14))
color_scheme_set("viridis")
dir.create("resultados", showWarnings = FALSE)

## ---- diagnostico-cargar ----
fit_data <- readRDS("data/ajuste_magallanes.rds")
chains <- fit_data$draws
config <- fit_data$configuration
stopifnot(length(dim(chains)) == 3, dim(chains)[2] == config$chains,
  dim(chains)[1] == config$iter - config$warmup, all(is.finite(chains)))
draws <- posterior::as_draws_array(chains)
parameters <- dimnames(chains)[[3]]
print(dim(chains)) # iteraciones posteriores x cadenas x variables

## ---- diagnostico-tabla ----
summary_table <- posterior::summarise_draws(draws,
  "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")
summary_table$revisar <- with(summary_table,
  !is.finite(rhat) | rhat >= 1.01 | ess_bulk < 400 | ess_tail < 400)
print(as.data.frame(summary_table), digits = 4)
# No retirar intercepto, efectos de comuna o escalas para mejorar el resumen.

## ---- diagnostico-nuts ----
nuts_table <- do.call(rbind, lapply(seq_along(fit_data$sampler), function(k) {
  s <- fit_data$sampler[[k]]
  energy <- s[, "energy__"]
  # Estimacion empirica de E-BFMI por cadena.
  bfmi <- sum(diff(energy)^2) / sum((energy - mean(energy))^2)
  data.frame(cadena = k, divergencias = sum(s[, "divergent__"]),
    limite_profundidad = sum(s[, "treedepth__"] >= config$max_treedepth),
    ebfmi = bfmi)
}))
nuts_table$revisar <- with(nuts_table,
  divergencias > 0 | limite_profundidad > 0 | !is.finite(ebfmi) | ebfmi < 0.3)
print(nuts_table, digits = 4)

## ---- diagnostico-trazas ----
# Incluir el parametro de peor R-hat junto a un coeficiente y una escala.
worst <- summary_table$variable[which.max(summary_table$rhat)]
selected <- unique(c(worst, "beta[1]", "sigma"))
g_trace <- mcmc_trace(draws, pars = selected,
  facet_args = list(ncol = 1)) +
  labs(colour = "Cadena", x = "Iteración posterior (sin calentamiento)",
       subtitle = "Magallanes: cadenas completas conservadas en su orden")
ggsave("resultados/trazas.png", g_trace, width = 12, height = 6, dpi = 120, bg = "white")

## ---- diagnostico-rangos ----
g_rank <- mcmc_rank_overlay(draws, pars = selected,
  facet_args = list(ncol = length(selected))) +
  labs(colour = "Cadena", subtitle = "Comparar cómo cada cadena ocupa los rangos de la misma variable")
ggsave("resultados/rangos.png", g_rank, width = 12, height = 4, dpi = 120, bg = "white")

## ---- diagnostico-acf ----
g_acf <- mcmc_acf(draws, pars = selected, lags = 30) +
  labs(subtitle = "Persistencia entre iteraciones: las simulaciones no son independientes")
ggsave("resultados/autocorrelacion.png", g_acf, width = 12, height = 5.5, dpi = 120, bg = "white")

## ---- diagnostico-pares ----
g_pair <- mcmc_scatter(draws, pars = c("beta[1]", "gamma[1]"), alpha = 0.12) +
  labs(x = "Intercepto", y = paste("Efecto de comuna", config$groups[1]),
       subtitle = "Asociacion posterior: revisar junto con R-hat, ESS y trazas")
ggsave("resultados/pares.png", g_pair, width = 8, height = 5, dpi = 120, bg = "white")

## ---- diagnostico-ejemplo ----
# Ejemplo artificial para aprender: no son las cadenas del modelo censal.
set.seed(2611)
example_draws <- array(NA_real_, c(500, 4, 2),
  dimnames = list(NULL, NULL, c("similares", "desplazada")))
example_draws[, , 1] <- matrix(rnorm(2000), 500, 4)
example_draws[, , 2] <- example_draws[, , 1]
example_draws[, 4, 2] <- example_draws[, 4, 2] + 2
example_summary <- posterior::summarise_draws(
  posterior::as_draws_array(example_draws), "rhat", "ess_bulk", "ess_tail")
g_example <- mcmc_trace(example_draws) +
  labs(colour = "Cadena", subtitle = "Ejemplo simulado: una cadena desplazada no debe ocultarse")
ggsave("resultados/ejemplo_simulado.png", g_example, width = 12, height = 5, dpi = 120, bg = "white")

## ---- diagnostico-replicas ----
# Replica MIXTA: nueva variacion de segmento, conservando X, Z y viviendas
# de los APTOS usados para ajustar. No es una validacion fuera de muestra.
chain_matrix <- posterior::as_draws_matrix(draws)
ids <- unique(round(seq(1, nrow(chain_matrix), length.out = min(500, nrow(chain_matrix)))))
beta_names <- paste0("beta[", seq_len(ncol(fit_data$X_obs)), "]")
gamma_names <- paste0("gamma[", seq_len(ncol(fit_data$Z_obs)), "]")
beta <- chain_matrix[ids, beta_names, drop = FALSE]
gamma <- chain_matrix[ids, gamma_names, drop = FALSE]
sigma <- as.numeric(chain_matrix[ids, "sigma"])
eta <- beta %*% t(fit_data$X_obs) + gamma %*% t(fit_data$Z_obs)
set.seed(2612)
density <- exp(eta + matrix(rnorm(length(eta)), nrow(eta)) * sigma)
# En el ajuste, la densidad NO tiene el techo predictivo delta_max.
# Se conservan los pisos/topes de lambda que si estan en la verosimilitud.
lambda <- pmin(pmax(sweep(density, 2, fit_data$observed$n_vpo, "*"), 1), 1e7)
y_rep <- matrix(rpois(length(lambda), lambda), nrow(lambda))
stopifnot(all(is.finite(y_rep)), all(y_rep >= 0))

## ---- diagnostico-chequeos ----
y_obs <- fit_data$observed$n_per
v_obs <- fit_data$observed$n_vpo
zones <- as.character(fit_data$observed$zona)
zone_check <- do.call(rbind, lapply(sort(unique(zones)), function(z) {
  columns <- zones == z
  rate <- rowSums(y_rep[, columns, drop = FALSE]) / sum(v_obs[columns])
  data.frame(zona = z, observado = sum(y_obs[columns]) / sum(v_obs[columns]),
    mediana_replica = median(rate), inferior = unname(quantile(rate, 0.025)),
    superior = unname(quantile(rate, 0.975)))
}))
print(zone_check)
g_ppc <- ggplot(zone_check, aes(y = zona)) +
  geom_segment(aes(x = inferior, xend = superior, yend = zona), colour = "#1f4e79", linewidth = 1) +
  geom_point(aes(x = mediana_replica), colour = "#1f4e79", size = 3) +
  geom_point(aes(x = observado), colour = "#c55a11", shape = 4, size = 4, stroke = 1.2) +
  labs(x = "Personas por vivienda: datos y replicas mixtas", y = NULL,
       subtitle = "Azul: mediana y rango 95 % simulado; cruz naranja: observado en aptos")
ggsave("resultados/chequeo_zonas.png", g_ppc, width = 11, height = 4.5, dpi = 120, bg = "white")
# Comparar tambien la dispersion de personas/vivienda entre segmentos.
observed_spread <- sd(y_obs / v_obs)
replicated_spread <- apply(sweep(y_rep, 2, v_obs, "/"), 1, sd)
spread_check <- data.frame(observado = observed_spread,
  inferior = unname(quantile(replicated_spread, 0.025)),
  mediana_replica = median(replicated_spread),
  superior = unname(quantile(replicated_spread, 0.975)))
print(spread_check)

## ---- diagnostico-guardar ----
write.csv(as.data.frame(summary_table), "resultados/diagnostico_parametros.csv", row.names = FALSE)
write.csv(nuts_table, "resultados/diagnostico_cadenas.csv", row.names = FALSE)
write.csv(zone_check, "resultados/chequeo_zonas.csv", row.names = FALSE)
saveRDS(list(parameters = as.data.frame(summary_table), nuts = nuts_table,
  example = as.data.frame(example_summary), zone_check = zone_check,
  spread_check = spread_check, worst = worst,
  sampler_review = any(summary_table$revisar) || any(nuts_table$revisar),
  configuration = config), "resultados/resumen.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")

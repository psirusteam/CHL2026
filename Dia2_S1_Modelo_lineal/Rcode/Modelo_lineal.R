##############################################################################
## Dia 2 - Sesion 1: Introduccion a modelos lineales con metodos bayesianos ##
## Regresion bayesiana de log(personas por vivienda) en Aysen.             ##
## Ejecutar con el directorio de trabajo en Dia2_S1_Modelo_lineal/:        ##
##   source("Rcode/Modelo_lineal.R")                                        ##
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

covariables <- grep("_std$", names(segmentos), value = TRUE)
covariables

## ---- construir-x ----
# Matriz de diseno: intercepto, zona (rural = referencia) y covariables estandarizadas
construir_X <- function(datos) {
  cbind(intercepto  = 1,
        zona_urbana = as.numeric(datos$zona == "urbana"),
        zona_aldea  = as.numeric(datos$zona == "aldea"),
        as.matrix(datos[, covariables]))
}
X <- construir_X(aptos)
stopifnot(all(is.finite(aptos$delta)), all(aptos$delta > 0),
          all(is.finite(X)), nrow(X) > ncol(X), qr(X)$rank == ncol(X))
y <- log(aptos$delta)
dim(X)

## ---- guardar-datos ----
guardar_valor("n_aptos", nrow(X))
guardar_valor("p", ncol(X))

## ---- exploracion ----
grafico_exploracion <- aptos %>%
  select(delta, zona, nlav_mean_std, dyn_world_urb_mean_std, ine_opb_area_ratio_std, dens_len_allroads_std) %>%
  pivot_longer(ends_with("_std"), names_to = "covariable", values_to = "valor") %>%
  mutate(covariable = str_remove(covariable, "_std$")) %>%
  ggplot(aes(valor, log(delta))) +
  geom_point(aes(colour = zona), alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black") +
  facet_wrap(~ covariable, scales = "free_x", nrow = 1) +
  labs(x = "Valor estandarizado", y = "log(personas por vivienda)", colour = "Zona") +
  theme(legend.position = "bottom")
guardar_figura(grafico_exploracion, "grafico_exploracion", alto = 5)

## ---- mco ----
ajuste_mco <- lm(y ~ X - 1)
tabla_mco <- broom::tidy(ajuste_mco) %>%
  mutate(term = str_remove(colnames(X), "_std$")) %>%
  select(termino = term, estimacion = estimate, error_est = std.error, valor_p = p.value)
tabla_mco

## ---- guardar-mco ----
guardar(tabla_mco, "tabla_mco")
guardar_valor("sigma_mco", sigma(ajuste_mco))

## ---- condicional-analitica ----
# Con beta ~ N(0, I) y sigma fija en su estimador, la condicional de beta es normal:
#   beta | sigma, y ~ N( (X'X/sigma^2 + I)^{-1} X'y/sigma^2 , (X'X/sigma^2 + I)^{-1} )
s <- sigma(ajuste_mco)
precision <- crossprod(X) / s^2 + diag(ncol(X))
covarianza <- solve(precision)
media_condicional <- covarianza %*% crossprod(X, y) / s^2

## ---- ajuste-stan ----
modelo_lineal <- stan_model("stan/regresion_lineal.stan")

datos_stan <- list(n = nrow(X), p = ncol(X), X = X, y = y, sd_beta = 1)

fit_lineal <- sampling(
  modelo_lineal, data = datos_stan,
  chains = 4, iter = 2000, warmup = 1000,
  seed = 1234, refresh = 0
)

## ---- diagnostico-rapido ----
resumen <- summarise_draws(as_draws_df(fit_lineal), "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
  filter(str_starts(variable, "beta|sigma")) %>%
  mutate(termino = c(colnames(X), "sigma")[match(variable,
    c(paste0("beta[", seq_len(ncol(X)), "]"), "sigma"))], .after = variable)

c(Rhat_max = max(resumen$rhat), ESS_min = min(resumen$ess_bulk))

## ---- guardar-diagnostico ----
guardar(tibble(Rhat_max = max(resumen$rhat), ESS_bulk_min = round(min(resumen$ess_bulk)),
               ESS_tail_min = round(min(resumen$ess_tail))), "diagnostico")

## ---- coeficientes ----
draws_beta <- as.data.frame(fit_lineal, pars = "beta") %>% setNames(colnames(X))

grafico_coeficientes <- mcmc_intervals(select(draws_beta, -intercepto), prob = 0.5, prob_outer = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Coeficiente (escala log de la densidad)")
guardar_figura(grafico_coeficientes, "grafico_coeficientes", alto = 6)

## ---- comparacion ----
comparacion <- resumen %>%
  filter(termino != "sigma") %>%
  transmute(termino = str_remove(termino, "_std$"),
            media_bayes = mean, de_bayes = sd,
            condicional_sigma_fija = media_condicional[, 1],
            mco = coef(ajuste_mco),
            efecto_pct = if_else(termino == "intercepto", NA_real_,
              colMeans(100 * (exp(as.matrix(draws_beta)) - 1))))
comparacion

## ---- guardar-comparacion ----
guardar(comparacion, "comparacion")

## ---- previa-vaga ----
fit_vaga <- sampling(modelo_lineal, data = modifyList(datos_stan, list(sd_beta = 100)),
                     chains = 4, iter = 2000, warmup = 1000, seed = 5678, refresh = 0)

grafico_previa <- tibble(termino = str_remove(colnames(X), "_std$"),
                         `N(0, 1)` = colMeans(draws_beta),
                         `N(0, 100²)` = colMeans(as.data.frame(fit_vaga, pars = "beta"))) %>%
  pivot_longer(-termino, names_to = "previa", values_to = "media") %>%
  ggplot(aes(media, fct_rev(factor(termino, levels = str_remove(colnames(X), "_std$"))), colour = previa)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(x = "Media posterior", y = NULL, colour = "Previa")
guardar_figura(grafico_previa, "grafico_previa", alto = 5.5)

## ---- r2 ----
sigma_draws <- as.data.frame(fit_lineal)$sigma
var_ajustada <- apply(as.matrix(draws_beta) %*% t(X), 1, var)
r2 <- var_ajustada / (var_ajustada + sigma_draws^2)
r2_mco <- 1 - sum(residuals(ajuste_mco)^2) / sum((y - mean(y))^2)   # lm sin intercepto no centra el R2

tabla_r2 <- tibble(R2_media = mean(r2), q2.5 = quantile(r2, 0.025), q97.5 = quantile(r2, 0.975),
                   R2_mco = r2_mco)
tabla_r2

## ---- guardar-r2 ----
guardar(tabla_r2, "tabla_r2")
grafico_r2 <- ggplot(tibble(r2), aes(r2)) +
  geom_histogram(bins = 40, fill = "#2e75b6", colour = "white") +
  geom_vline(xintercept = r2_mco, colour = "#c00000", linewidth = 1.2) +
  labs(x = expression(R^2~bayesiano), y = "Draws", subtitle = "Línea roja: R² de mínimos cuadrados")
guardar_figura(grafico_r2, "grafico_r2", alto = 4)

## ---- ppc ----
y_rep <- as.matrix(fit_lineal, pars = "y_rep")
set.seed(1)
idx <- sample(nrow(y_rep), 50)

grafico_ppc <- ppc_dens_overlay(y, y_rep[idx, ]) + labs(x = "log(personas por vivienda)") +
  ppc_stat_grouped(y, y_rep, group = aptos$zona, stat = "sd") +
  plot_layout(widths = c(1, 1.3))
guardar_figura(grafico_ppc, "grafico_ppc", alto = 4.8)

## ---- heterocedasticidad ----
grafico_residuos <- aptos %>%
  mutate(residuo = y - fitted(ajuste_mco),
         viviendas = cut(n_vpo, c(0, 5, 10, 20, 40, Inf),
                         labels = c("1–5", "6–10", "11–20", "21–40", "más de 40"))) %>%
  ggplot(aes(viviendas, residuo)) +
  geom_boxplot(fill = "#9ecae1", outlier.alpha = 0.3) +
  labs(x = "Viviendas ocupadas del segmento", y = "Residuo de log(δ)")
tabla_residuos <- aptos %>%
  mutate(residuo = y - fitted(ajuste_mco),
         viviendas = cut(n_vpo, c(0, 5, 10, 20, 40, Inf),
                         labels = c("1–5", "6–10", "11–20", "21–40", "más de 40"))) %>%
  group_by(viviendas) %>%
  summarise(segmentos = n(), desv_residuo = sd(residuo),
            desv_poisson_aprox = mean(sqrt(1 / (n_vpo * delta))))

## ---- guardar-heterocedasticidad ----
guardar_figura(grafico_residuos, "grafico_residuos", alto = 4.5)
guardar(tabla_residuos, "tabla_residuos")

## ---- prediccion ----
X_todos <- construir_X(segmentos)
stopifnot(identical(colnames(X_todos), colnames(X)), all(is.finite(X_todos)))
set.seed(1)
eta <- as.matrix(draws_beta) %*% t(X_todos)                        # draws x segmentos
delta_pred <- exp(eta + matrix(rnorm(length(eta)), nrow(eta)) * sigma_draws)

tabla_prediccion <- segmentos %>%
  mutate(delta_pred = colMeans(delta_pred)) %>%
  group_by(clasificacion_cobertura) %>%
  summarise(segmentos = n(), delta_observado = mean(delta[is.finite(delta)], na.rm = TRUE),
            delta_predicho = mean(delta_pred))
tabla_prediccion

## ---- guardar-prediccion ----
guardar(tabla_prediccion, "tabla_prediccion")

## ---- diagnostics-both-fits ----
fit_list <- list(regularized = fit_lineal, wide_prior = fit_vaga)
diagnostics <- purrr::imap_dfr(fit_list, function(fit, name) {
  posterior::summarise_draws(posterior::as_draws_array(fit),
                           "rhat", "ess_bulk", "ess_tail", "mcse_mean") %>%
    filter(grepl("^(beta\\[|sigma$)", variable)) %>%
    mutate(fit = name, .before = 1)
})
guardar(diagnostics, "diagnostico_completo")
nuts_diagnostics <- purrr::imap_dfr(fit_list, function(fit, name) {
  chains <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  tibble(fit = name,
         divergences = sum(vapply(chains, function(x) sum(x[, "divergent__"]), numeric(1))),
         max_treedepth_hits = sum(vapply(chains, function(x) sum(x[, "treedepth__"] >= 10), numeric(1))))
})
guardar(nuts_diagnostics, "diagnostico_nuts")
writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")


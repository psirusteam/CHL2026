# D3S3: modelo de poblacion de Chile, demostracion guiada en Aysen.
# Ejecutar por bloques desde la carpeta de esta sesion.
# Copia exacta del Stan del proyecto a la fecha de revision; no modificar
# el archivo nacional desde la carpeta del taller.

## ---- chile-paquetes ----
library(rstan)
library(posterior)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))
dir.create("resultados", showWarnings = FALSE)
theme_set(theme_minimal(base_size = 15))

## ---- chile-datos ----
segmentos <- readRDS("data/segmentos_aysen.rds")
stopifnot(!anyDuplicated(segmentos[c("edid", "fid")]),
          !anyNA(segmentos[c("zona", "cut", "n_vpo", "n_per", "apto_modelo", "excluido_colectivo")]),
          all(is.finite(segmentos$n_vpo)), all(segmentos$n_vpo >= 0),
          all(segmentos$n_vpo == floor(segmentos$n_vpo)),
          all(is.finite(segmentos$n_per)), all(segmentos$n_per >= 0),
          all(segmentos$n_per == floor(segmentos$n_per)))
# Respetar la elegibilidad preparada y verificar exposicion utilizable.
aptos <- segmentos[segmentos$apto_modelo, ]
stopifnot(all(aptos$n_vpo > 0), all(!aptos$excluido_colectivo))
# Prediccion didactica: viviendas particulares con VPO positiva.
# No es el universo completo regional ni un total poblacional validado.
destino <- segmentos[!segmentos$excluido_colectivo & segmentos$n_vpo > 0, ]
resumen_datos <- data.frame(
  conjunto = c("Base disponible", "Aptos para ajustar", "Destino de la demostracion"),
  segmentos = c(nrow(segmentos), nrow(aptos), nrow(destino)))
print(resumen_datos)

## ---- chile-matrices ----
# Seleccion pequena para aprender el flujo: no reproduce el ajuste nacional.
# Los _std ya usan los parametros nacionales; no reestandarizar por region.
covariables <- c("nlav_mean_std", "dyn_world_urb_mean_std")
construir_X <- function(d) {
  cbind(intercepto = 1, zona_aldea = as.numeric(d$zona == "aldea"),
        zona_urbana = as.numeric(d$zona == "urbana"),
        as.matrix(d[, covariables, drop = FALSE]))
}
# El diccionario se define UNA vez para conservar columnas y orden.
niveles_comuna <- sort(unique(as.character(segmentos$cut)))
construir_Z <- function(d) {
  z <- outer(as.character(d$cut), niveles_comuna, "==") * 1
  colnames(z) <- niveles_comuna
  z
}
X_obs <- construir_X(aptos)
X_tot <- construir_X(destino)
Z_obs <- construir_Z(aptos)
Z_tot <- construir_Z(destino)
stopifnot(all(is.finite(X_obs)), all(is.finite(X_tot)),
          qr(X_obs)$rank == ncol(X_obs),
          identical(colnames(X_obs), colnames(X_tot)),
          identical(colnames(Z_obs), colnames(Z_tot)),
          all(rowSums(Z_obs) == 1), all(rowSums(Z_tot) == 1))
print(data.frame(comuna = niveles_comuna, aptos = colSums(Z_obs)))

## ---- chile-lista ----
# Escenario base: usar VPO observada en prediccion, sin ampliacion de cobertura.
cobertura_destino <- rep(1, nrow(destino))
# Valor de referencia historico del taller (p99 nacional de los aptos).
# No se estima en Aysen ni representa un limite universal biologico.
delta_max <- 4.35
datos_stan <- list(D_obs = nrow(aptos), D_tot = nrow(destino),
  P = ncol(X_obs), Q = ncol(Z_obs),
  Y_obs = as.integer(aptos$n_per), V_obs = as.integer(aptos$n_vpo),
  V_tot = as.integer(round(destino$n_vpo / cobertura_destino)),
  X_obs = X_obs, X_tot = X_tot, Z_obs = Z_obs, Z_tot = Z_tot,
  delta_max = delta_max)
str(datos_stan[c("D_obs", "D_tot", "P", "Q", "delta_max")])

## ---- chile-ajuste ----
modelo <- stan_model("stan/modelo_poblacion_v02.stan")
# Valores iniciales razonables: no son distribuciones previas.
inicio <- lm(log(n_per / n_vpo) ~ X_obs - 1, data = aptos)
inicializar <- function() list(
  beta = unname(coef(inicio)) + rnorm(ncol(X_obs), 0, 0.01),
  gamma = rep(0, ncol(Z_obs)), sigmaU = 0.1,
  sigma = max(sd(residuals(inicio)), 0.1),
  densidad = aptos$n_per / aptos$n_vpo * exp(rnorm(nrow(aptos), 0, 0.01)))
ajuste <- sampling(modelo, data = datos_stan, init = inicializar,
  chains = 4, iter = 3000, warmup = 1500, seed = 2533, refresh = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  pars = c("beta", "gamma", "sigma", "sigmaU", "Y_tot"))

## ---- chile-diagnostico ----
# Revisar TODOS los coeficientes y escalas: no excluir el intercepto/gamma.
diagnostico <- posterior::summarise_draws(
  posterior::as_draws_array(as.array(ajuste, pars = c("beta", "gamma", "sigma", "sigmaU"))),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail")
sampler <- get_sampler_params(ajuste, inc_warmup = FALSE)
controles <- data.frame(
  rhat_max = max(diagnostico$rhat),
  ess_bulk_min = min(diagnostico$ess_bulk), ess_tail_min = min(diagnostico$ess_tail),
  divergencias = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
  limite_profundidad = sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 12), numeric(1))))
controles$revisar <- with(controles, !is.finite(rhat_max) | rhat_max >= 1.01 |
  ess_bulk_min < 400 | ess_tail_min < 400 | divergencias > 0 | limite_profundidad > 0)
print(controles, digits = 5)
if (controles$revisar) warning("El ajuste requiere revision: no utilizar sus predicciones para decisiones.")

## ---- chile-resultados ----
beta <- as.matrix(ajuste, pars = "beta")
factores <- data.frame(termino = colnames(X_obs)[-1],
  t(apply(exp(beta[, -1, drop = FALSE]), 2, quantile, c(0.025, 0.5, 0.975))),
  check.names = FALSE)
print(factores)
# Un coeficiente de covariable: factor multiplicativo por una DE nacional,
# manteniendo zona, comuna y restantes covariables. No implica causalidad.
y_tot <- as.matrix(ajuste, pars = "Y_tot")
stopifnot(all(is.finite(y_tot)), all(y_tot >= 0))
prediccion <- data.frame(destino[c("edid", "fid", "cut", "zona", "apto_modelo", "n_per", "n_vpo")],
  media = colMeans(y_tot),
  t(apply(y_tot, 2, quantile, c(0.025, 0.5, 0.975))), check.names = FALSE)
# Guardar por draw para sumar antes de resumir: uso en sesiones posteriores.
total_draw <- rowSums(y_tot)
resumen_total <- data.frame(universo = "Destino didactico: no es toda la region",
  censado_mismo_subconjunto = sum(destino$n_per),
  media_predictiva = mean(total_draw),
  inferior = unname(quantile(total_draw, 0.025)),
  superior = unname(quantile(total_draw, 0.975)), revisar_ajuste = controles$revisar)
print(resumen_total)

## ---- chile-figuras ----
g_coef <- ggplot(factores, aes(x = `50%`, y = termino)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, yend = termino)) +
  geom_point(size = 2.8, colour = "#1f4e79") +
  scale_y_discrete(labels = c(zona_aldea = "Aldea frente a rural", zona_urbana = "Urbana frente a rural",
    nlav_mean_std = "Luces nocturnas", dyn_world_urb_mean_std = "Cobertura urbana")) +
  labs(x = "Factor multiplicativo: mediana e intervalo del 95 %", y = NULL,
       subtitle = "Aysén: comparación por zona; covariables por una desviación estándar nacional")
ggplot2::ggsave("resultados/factores.png", g_coef, width = 11, height = 4, dpi = 120, bg = "white")

## ---- chile-guardar ----
write.csv(resumen_datos, "resultados/resumen_datos.csv", row.names = FALSE)
write.csv(as.data.frame(diagnostico), "resultados/diagnostico.csv", row.names = FALSE)
write.csv(prediccion, "resultados/prediccion_segmentos.csv", row.names = FALSE)
saveRDS(list(resumen_datos = resumen_datos, controles = controles,
  factores = factores, resumen_total = resumen_total,
  configuracion = list(covariables = covariables, columnas_X = colnames(X_obs),
    niveles_comuna = niveles_comuna, cobertura = 1, delta_max = delta_max,
    semilla = 2533, cadenas = 4, iter = 3000, warmup = 1500,
    firma_stan = unname(tools::md5sum("stan/modelo_poblacion_v02.stan")))),
  "resultados/resumen.rds")
saveRDS(total_draw, "resultados/total_draw_destino.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")
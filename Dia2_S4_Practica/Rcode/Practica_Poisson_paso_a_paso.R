# Practica D2S4. Ejecutar bloque por bloque desde la carpeta de la sesion.
# Este script ajusta el modelo Poisson por zona y contiene todos los pasos de la sesion.

## ---- taller-paquetes ----
library(rstan)
library(posterior)
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))
dir.create("resultados/practica_poisson", recursive = TRUE, showWarnings = FALSE)

## ---- taller-datos ----
censo <- readRDS("data/segmentos_magallanes.rds")
table(censo$zona, censo$apto_modelo)
# La bandera de aptitud viene de la preparacion de los datos del taller.
# Revisar sus criterios con el docente; no equivale a representar todo el censo.
datos <- censo[!is.na(censo$apto_modelo) & censo$apto_modelo, ]
datos$zona <- factor(datos$zona, levels = c("rural", "aldea", "urbana"))
stopifnot(!anyNA(datos[c("n_per", "n_vpo", "zona")]),
          all(is.finite(datos$n_per)), all(datos$n_per >= 0),
          all(datos$n_per == floor(datos$n_per)),
          all(is.finite(datos$n_vpo)), all(datos$n_vpo > 0))
aggregate(cbind(n_per, n_vpo) ~ zona, datos, sum)

## ---- taller-lista ----
# Rural es la referencia. Las otras columnas indican aldea y urbana.
X <- model.matrix(~ zona, data = datos)
head(X)
datos_stan <- list(n = nrow(datos), p = ncol(X), X = X,
                   y = as.integer(datos$n_per), log_v = log(datos$n_vpo))
stopifnot(qr(X)$rank == ncol(X))

## ---- taller-ajuste ----
# Compilar transforma el archivo Stan en un programa ejecutable.
modelo <- stan_model("stan/poisson_offset.stan")
# Cuatro cadenas: cuatro recorridos para explorar los parametros posibles.
# Las primeras 1000 iteraciones de cada cadena son de adaptacion.
beta_inicial <- c(log(sum(datos$n_per) / sum(datos$n_vpo)), rep(0, ncol(X) - 1))
ajuste <- sampling(modelo, data = datos_stan,
  init = function() list(beta = beta_inicial + rnorm(ncol(X), 0, 0.05)),
  chains = 4, iter = 2000, warmup = 1000, seed = 2411, refresh = 500,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  pars = c("beta", "y_rep", "rng_out_of_range"))

## ---- taller-diagnostico ----
if (!interactive()) pdf("resultados/practica_poisson/diagnosticos.pdf", width = 12, height = 7)
diagnostico <- posterior::summarise_draws(
  posterior::as_draws_array(as.array(ajuste, pars = "beta")),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail")
print(as.data.frame(diagnostico), digits = 4)
sampler <- get_sampler_params(ajuste, inc_warmup = FALSE)
divergencias <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
limite_profundidad <- sum(vapply(sampler,
  function(x) sum(x[, "treedepth__"] >= 12), numeric(1)))
print(c(divergencias = divergencias, limite_profundidad = limite_profundidad))
traceplot(ajuste, pars = "beta")
# No interpretar si fallan estos controles: pedir ayuda al docente.
stopifnot(all(is.finite(diagnostico$rhat)), all(diagnostico$rhat < 1.01),
          all(diagnostico$ess_bulk >= 400), all(diagnostico$ess_tail >= 400),
          divergencias == 0, limite_profundidad == 0,
          all(as.matrix(ajuste, pars = "rng_out_of_range") == 0))

## ---- taller-interpretacion ----
beta <- as.matrix(ajuste, pars = "beta")
# Transformar coeficientes a personas esperadas por vivienda, por zona.
X_zona <- model.matrix(~ zona, data = data.frame(
  zona = factor(c("rural", "aldea", "urbana"), levels = levels(datos$zona))))
densidad <- exp(beta %*% t(X_zona))
resumen_zona <- data.frame(zona = c("rural", "aldea", "urbana"),
  t(apply(densidad, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  check.names = FALSE)
print(resumen_zona)
# Comparacion descriptiva con los datos usados para ajustar, no validacion externa.
observado_zona <- aggregate(cbind(n_per, n_vpo) ~ zona, datos, sum)
observado_zona$personas_por_vivienda <- observado_zona$n_per / observado_zona$n_vpo
print(observado_zona)

## ---- taller-prediccion ----
# Viviendas hipoteticas de zona urbana: cambiar 100 por 150.
viviendas_nuevas <- 100
set.seed(2412)
media_esperada <- viviendas_nuevas * densidad[, 3]
poblacion_simulada <- rpois(length(media_esperada), media_esperada)
print(rbind(
  media_esperada = quantile(media_esperada, c(0.025, 0.5, 0.975)),
  conteo_predictivo = quantile(poblacion_simulada, c(0.025, 0.5, 0.975))))

## ---- taller-chequeo ----
# Un conteo replicado por segmento y simulacion bajo el modelo ajustado.
y_rep <- as.matrix(ajuste, pars = "y_rep")
stopifnot(all(y_rep >= 0))
mu <- colMeans(exp(beta %*% t(X))) * datos$n_vpo
# Indicador preparado: comparar variacion observada con variacion simulada.
indicador <- function(y) sum((y - mu)^2 / mu) / (length(y) - ncol(X))
observado <- indicador(datos$n_per)
ids <- unique(round(seq(1, nrow(y_rep), length.out = min(500, nrow(y_rep)))))
simulado <- apply(y_rep[ids, , drop = FALSE], 1, indicador)
print(c(observado = observado, quantile(simulado, c(0.05, 0.5, 0.95))))
hist(simulado, breaks = 30, col = "lightblue",
  main = "Variacion reproducida por el Poisson", xlab = "Indicador de variacion",
  xlim = range(c(simulado, observado)))
abline(v = observado, col = "red", lwd = 2)
legend("topright", "Datos observados", col = "red", lwd = 2, bty = "n")

## ---- taller-guardar ----
write.csv(resumen_zona, "resultados/practica_poisson/densidad_zona.csv", row.names = FALSE)
saveRDS(list(diagnostico = diagnostico, divergencias = divergencias,
  limite_profundidad = limite_profundidad, densidad_zona = resumen_zona,
  viviendas_nuevas = viviendas_nuevas,
  prediccion = quantile(poblacion_simulada, c(0.025, 0.5, 0.975)),
  dispersion_observada = observado,
  dispersion_simulada = quantile(simulado, c(0.05, 0.5, 0.95)), seed = 2411),
  "resultados/practica_poisson/resumen.rds")
capture.output(sessionInfo(), file = "resultados/practica_poisson/sessionInfo.txt")
if (!interactive()) dev.off()

## ---- taller-pln-ajuste ----
modelo_pln <- stan_model("stan/poisson_lognormal.stan")
beta_pln_inicial <- colMeans(beta)
ajuste_pln <- sampling(modelo_pln, data = datos_stan,
  init = function() list(beta = beta_pln_inicial + rnorm(ncol(X), 0, 0.02),
                        sigma = 0.15, z = rep(0, nrow(X))),
  chains = 4, iter = 2000, warmup = 1000, seed = 2413, refresh = 500,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  pars = c("beta", "sigma", "z", "y_rep", "rng_out_of_range"))

## ---- taller-pln-diagnostico ----
parametros_pln <- as.array(ajuste_pln, pars = c("beta", "sigma", "z"))
diagnostico_pln <- posterior::summarise_draws(
  posterior::as_draws_array(parametros_pln), "mean", "sd", "rhat", "ess_bulk", "ess_tail")
sampler_pln <- get_sampler_params(ajuste_pln, inc_warmup = FALSE)
control_pln <- c(Rhat_max = max(diagnostico_pln$rhat),
  ESS_bulk_min = min(diagnostico_pln$ess_bulk),
  ESS_tail_min = min(diagnostico_pln$ess_tail),
  divergencias = sum(vapply(sampler_pln, function(x) sum(x[, "divergent__"]), numeric(1))),
  limite_profundidad = sum(vapply(sampler_pln,
    function(x) sum(x[, "treedepth__"] >= 12), numeric(1))))
print(round(control_pln, 3))
print(as.data.frame(diagnostico_pln[diagnostico_pln$variable == "sigma", ]), digits = 3)
stopifnot(all(is.finite(control_pln)), control_pln["Rhat_max"] < 1.01,
          control_pln["ESS_bulk_min"] >= 400, control_pln["ESS_tail_min"] >= 400,
          control_pln["divergencias"] == 0, control_pln["limite_profundidad"] == 0,
          all(as.matrix(ajuste_pln, pars = "rng_out_of_range") == 0))

## ---- taller-pln-chequeo ----
y_rep_pln <- as.matrix(ajuste_pln, pars = "y_rep")
stopifnot(all(y_rep_pln >= 0))
# El mismo indicador, los mismos segmentos y la misma referencia Poisson.
simulado_pln <- apply(y_rep_pln[ids, , drop = FALSE], 1, indicador)
intervalos <- rbind(Poisson = quantile(simulado, c(0.05, 0.5, 0.95)),
                    `Poisson–lognormal` = quantile(simulado_pln, c(0.05, 0.5, 0.95)))
print(intervalos)
print(observado)

## ---- taller-pln-guardar ----
ruta_pln <- "resultados/practica_poisson"
graficar_comparacion <- function() {
  limites <- range(simulado, simulado_pln, observado)
  cortes <- seq(limites[1], limites[2], length.out = 33)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (modelo in rownames(intervalos)) {
    valores <- if (modelo == "Poisson") simulado else simulado_pln
    hist(valores, breaks = cortes, xlim = limites, col = "lightblue",
         main = modelo, xlab = "Indicador de variación", ylab = "Réplicas")
    abline(v = observado, col = "firebrick", lwd = 3)
  }
}
png(file.path(ruta_pln, "comparacion_pln.png"), width = 1380, height = 520, res = 120)
graficar_comparacion()
invisible(dev.off())
if (interactive()) graficar_comparacion()
saveRDS(list(intervalos = intervalos, observado = observado, control = control_pln,
             sigma = quantile(as.matrix(ajuste_pln, pars = "sigma"), c(0.05, 0.5, 0.95))),
        file.path(ruta_pln, "comparacion_pln.rds"))

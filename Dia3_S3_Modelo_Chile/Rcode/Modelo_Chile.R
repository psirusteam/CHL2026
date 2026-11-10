# Modelo final implementado en Chile: Aysén, 40 covariables.
# Ejecutar por bloques desde esta sesion. No modifica los ajustes fuente.

## ---- chile-paquetes ----
library(posterior)
library(ggplot2)
source("../utilidades_modelo_final.R")
dir.create("resultados", showWarnings = FALSE)
theme_set(theme_minimal(base_size = 15))

## ---- chile-datos ----
final <- load_final_model("11")
segmentos <- final$segments
aptos <- final$observed
destino <- segmentos # Mismas llaves y universo del ajuste final.
resumen_datos <- data.frame(conjunto = c("Base disponible", "Aptos para ajustar", "Destino del modelo final"),
 segmentos = c(nrow(segmentos), nrow(aptos), nrow(destino)))
print(resumen_datos)

## ---- chile-matrices ----
covariables <- final$covariates
# Construir el mismo diseno final y cotejarlo con los insumos guardados.
datos_stan <- final$data_stan
niveles_comuna <- final$configuration$group_labels
construir_X <- function(d) cbind(intercepto = 1,
 zona_urbana = as.numeric(d$zona == "urbana"),
 zona_aldea = as.numeric(d$zona == "aldea"),
 as.matrix(d[, covariables, drop = FALSE]))
construir_Z <- function(d) {
 z <- outer(as.character(d$cut), niveles_comuna, "==") * 1
 colnames(z) <- paste0("cut_", niveles_comuna)
 z
}
X_obs <- construir_X(aptos); X_tot <- construir_X(destino)
Z_obs <- construir_Z(aptos); Z_tot <- construir_Z(destino)
stopifnot(isTRUE(all.equal(X_obs, datos_stan$X_obs)),
 isTRUE(all.equal(X_tot, datos_stan$X_tot)),
 isTRUE(all.equal(Z_obs, datos_stan$Z_obs)),
 isTRUE(all.equal(Z_tot, datos_stan$Z_tot)))
stopifnot(ncol(X_obs) == 43L, length(covariables) == 40L,
 identical(colnames(X_obs), colnames(X_tot)),
 identical(colnames(Z_obs), colnames(Z_tot)), all(rowSums(Z_obs) == 1),
 all(rowSums(Z_tot) == 1), all(is.finite(X_tot)))
print(data.frame(comuna = niveles_comuna, aptos = colSums(Z_obs)))
# El rango y la colinealidad se inspeccionan, sin retirar covariables.
print(c(columnas_X = ncol(X_obs), rango_X = qr(X_obs)$rank))

## ---- chile-lista ----
delta_max <- datos_stan$delta_max
# V_tot ya contiene la regla de exposicion implementada en el proyecto.
# No reemplazar por n_vpo ni asumir cobertura uno como escenario final.
stopifnot(identical(as.numeric(destino$V_tot), as.numeric(datos_stan$V_tot)))
str(datos_stan[c("D_obs", "D_tot", "P", "Q", "delta_max")])
print(final$original_configuration)

## ---- chile-ajuste ----
# La clase usa las cadenas originales del ajuste implementado.
# Para repetir el muestreo, ejecutar este bloque con reestimar = TRUE.
reestimar <- FALSE
if (reestimar) {
 library(rstan)
 ajuste <- rstan::stan(file = "stan/modelo_poblacion_v02.stan", data = datos_stan,
  chains = final$configuration$chains, iter = final$configuration$iter,
  warmup = final$configuration$warmup, seed = final$configuration$seed_Stan,
  control = list(max_treedepth = final$configuration$max_treedepth))
 # Una repeticion es una corrida nueva: conservarla separada de la referencia.
 saveRDS(ajuste, "resultados/ajuste_repetido.rds")
}
# El resto del material analiza la referencia final, no la repeticion opcional.
chain_matrix <- posterior::as_draws_matrix(posterior::as_draws_array(final$draws))

## ---- chile-diagnostico ----
checks <- check_final_model(final)
diagnostico <- checks$diagnostics
controles <- checks$controls
print(controles, digits = 5)
#if (controles$revisar) warning("AJUSTE FINAL PENDIENTE DE VALIDACION: resultados solo para ejercicio de calculo.")

## ---- chile-resultados ----
beta <- chain_matrix[, paste0("beta[", seq_len(ncol(X_obs)), "]"), drop = FALSE]
factores <- data.frame(termino = colnames(X_obs)[-1],
  t(apply(exp(beta[, -1, drop = FALSE]), 2, quantile, c(0.025, 0.5, 0.975))),
  check.names = FALSE)
print(factores)
# Un coeficiente de covariable: factor multiplicativo por una DE nacional,
# manteniendo zona, comuna y restantes covariables. No implica causalidad.
y_tot <- final$Y_tot # Cantidades generadas originales, mismo orden de draws.
stopifnot(all(is.finite(y_tot)), all(y_tot >= 0))
prediccion <- data.frame(estado = "Pendiente de validacion", destino[c("edid", "fid", "cut", "zona", "apto_modelo", "n_per", "n_vpo")],
  media = colMeans(y_tot),
  t(apply(y_tot, 2, quantile, c(0.025, 0.5, 0.975))), check.names = FALSE)
# Guardar por draw para sumar antes de resumir: uso en sesiones posteriores.
total_draw <- rowSums(y_tot)
resumen_total <- data.frame(universo = "Universo regional del modelo final; pendiente de validacion",
  censado_mismo_subconjunto = sum(destino$n_per),
  media_predictiva = mean(total_draw),
  inferior = unname(quantile(total_draw, 0.025)),
  superior = unname(quantile(total_draw, 0.975)), revisar_ajuste = controles$revisar)
resumen_total$estado <- "Pendiente de validacion del muestreo"
print(resumen_total)

## ---- chile-figuras ----
write.csv(factores, "resultados/factores_completos.csv", row.names = FALSE)
# Mostrar algunos terminos; el ajuste conserva los 43 coeficientes.
plot_factors <- factores[factores$termino %in%
                           c("zona_urbana", "zona_aldea",
                             "nlav_mean_std",
                             "dyn_world_urb_mean_std", 
                             "ine_opb_mean_barea_std", 
                             "mean_slope_std"), ]
g_coef <- ggplot(plot_factors, aes(x = `50%`, y = termino)) +
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
    niveles_comuna = niveles_comuna, exposicion = "V_tot del modelo final", delta_max = delta_max,
    semilla = final$configuration$seed_Stan, cadenas = final$configuration$chains,
    iter = final$configuration$iter, warmup = final$configuration$warmup,
    firma_stan = unname(tools::md5sum("stan/modelo_poblacion_v02.stan")))),
  "resultados/resumen.rds")
saveRDS(total_draw, "resultados/total_draw_destino.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")

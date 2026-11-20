# D4S2 | Predicción por segmento y por comuna en Magallanes
#
# Punto de partida: el ajuste final que revisamos en D4S1 y los segmentos
# que el proyecto definió como destino. Ejecutar por bloques desde
# Dia4_S2_Prediccion, manteniendo D4S1 en la carpeta hermana. Usamos las
# cadenas guardadas y generamos conteos nuevos en R.
#
# Recorrido: verificar insumos → simular personas por vivienda → obtener
# conteos para tres supuestos de cobertura → sumar por comuna y escenario.
# Pregunta de trabajo: ¿cuánto cambia el total del mismo destino cuando
# cambiamos las viviendas usadas como exposición?
#
# Para probar otra hipótesis, editar prediccion-escenarios. Para revisar otra
# comuna, cambiar selected_commune en prediccion-agregar. Para trabajar con
# otra región, usar Predecir_una_region.R con su paquete de insumos.

## ---- prediccion-paquetes ----
# posterior convierte las cadenas en una matriz; ggplot2 prepara las figuras.
library(posterior)
library(ggplot2)
theme_set(theme_minimal(base_size = 14))
# Ubicar la sesión al ejecutar desde el proyecto, la sesión o su carpeta Rcode.
candidatas <- c(".", "Dia4_S2_Prediccion", "..")
sesion <- candidatas[file.exists(file.path(candidatas, "data", "destino_magallanes.rds"))]
if (!length(sesion)) {
  stop("Abra el proyecto CHL2026ModelosPob o fije la carpeta Dia4_S2_Prediccion como directorio de trabajo.")
}
setwd(sesion[1])
source("Rcode/funciones_prediccion.R", local = TRUE, encoding = "UTF-8")
dir.create("resultados", showWarnings = FALSE)

## ---- prediccion-insumos ----
# fit_data trae parámetros ajustados en D4S1; target trae el destino a predecir.
# X_tot contiene zona y covariables; Z_tot indica la comuna de cada segmento.
# La firma del modelo, las llaves y el orden de las columnas deben coincidir.

fit_data <- readRDS("../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")
diagnostics <- readRDS("../Dia4_S1_Diagnostico/resultados/resumen.rds")
target <- readRDS("data/destino_magallanes.rds")
destino <- target$destino
X_tot <- target$X_tot
Z_tot <- target$Z_tot
config <- fit_data$configuration
stopifnot(identical(config, diagnostics$configuration),
  identical(target$stan_hash, config$stan_hash),
  identical(target$source_configuration, config),
  identical(colnames(X_tot), config$columns_X),
  identical(colnames(Z_tot), config$groups),
  !anyDuplicated(destino[c("edid", "fid")]),
  nrow(destino) == nrow(X_tot), nrow(destino) == nrow(Z_tot),
  all(is.finite(X_tot)), all(is.finite(Z_tot)),
  all(rowSums(Z_tot) == 1), all(Z_tot %in% c(0,1)),
  all(as.character(destino$cut) %in% config$group_labels),
  identical(as.integer(max.col(Z_tot, ties.method = "first")),
    as.integer(match(as.character(destino$cut), config$group_labels))),
  all(destino$n_vpo >= 0), all(is.finite(destino$V_tot)), all(destino$V_tot >= 0))

# Encontrar cada segmento apto en el destino y verificar que sus datos coincidan.

idx <- match(key(fit_data$observed), key(destino))
stopifnot(!anyNA(idx),
  identical(as.numeric(destino$n_per[idx]), as.numeric(fit_data$observed$n_per)),
  identical(as.numeric(destino$n_vpo[idx]), as.numeric(fit_data$observed$n_vpo)),
  isTRUE(all.equal(unname(X_tot[idx, ]), unname(fit_data$X_obs))),
  isTRUE(all.equal(unname(Z_tot[idx, ]), unname(fit_data$Z_obs))))

cat("Destino: mismo universo y V_tot del modelo final. Resultados pendientes de validacion.\n")

## ---- prediccion-draws ----
# Cada fila reúne beta, los efectos de comuna (gamma) y la dispersión sigma.
# Usar esa misma fila en todos los segmentos de una réplica conserva lo compartido.
# Tomar hasta 1.000 filas espaciadas para limitar memoria al predecir.
# R-hat y ESS se calcularon con las cadenas completas en D4S1.

all_draws <- posterior::as_draws_matrix(posterior::as_draws_array(fit_data$draws))

ids <- unique(round(seq(1, nrow(all_draws), length.out = min(1000, nrow(all_draws)))))
beta_names <- paste0("beta[", seq_len(ncol(X_tot)), "]")
gamma_names <- paste0("gamma[", seq_len(ncol(Z_tot)), "]")
beta <- all_draws[ids, beta_names, drop = FALSE]
gamma <- all_draws[ids, gamma_names, drop = FALSE]
sigma <- as.numeric(all_draws[ids, "sigma"])
eta <- beta %*% t(X_tot) + gamma %*% t(Z_tot)

set.seed(2621)
density_raw <- exp(eta + matrix(rnorm(length(eta)), nrow(eta)) * sigma)
delta_max <- config$delta_max
density <- pmin(density_raw, delta_max)
# Reutilizar los mismos uniformes para comparar las exposiciones bajo la misma simulación.
set.seed(2622)
uniforms <- matrix(runif(length(density)), nrow(density))

## ---- prediccion-escenarios ----
# Comparar el mismo destino: base, 95 % en todos y 90 % solo en Incompleta.
# Son valores de cobertura fijados para el ejercicio.

c_base <- rep(1, nrow(destino))
c_uniforme <- rep(0.95, nrow(destino))
c_focalizada <- ifelse(as.character(destino$clasificacion_cobertura) == "Incompleta", 0.90, 1)
coverage_scenarios <- list(base = c_base, uniforme_095 = c_uniforme,
                           focalizada_090 = c_focalizada)

## ---- prediccion-llaves ----
# La tabla externa podría llegar desordenada. Ordenarla según edid y fid del destino.
# La comprobación final debe recuperar exactamente el escenario focalizado.
# Ejemplo de archivo externo desordenado: alinear por edid/fid, nunca por posicion.

set.seed(2623)
coverage_table <- data.frame(destino[c("edid", "fid")], cobertura = c_focalizada)
coverage_table <- coverage_table[sample(nrow(coverage_table)), ]
stopifnot(identical(as.numeric(alinear_cobertura(coverage_table, destino)), as.numeric(c_focalizada)))

## ---- prediccion-resumenes ----
# Por segmento se resumen intensidad y conteo por separado.
# Por comuna y para el destino completo, primero sumar conteos en cada réplica.
# Comparar el conteo censal usando las mismas filas que entraron en la predicción.

scenario_tables <- segment_tables <- commune_tables <- list()
total_draws <- commune_draws <- list()
base_y <- base_lambda <- NULL

for (name in names(coverage_scenarios)) {
  pred <- predecir(coverage_scenarios[[name]], destino, density, uniforms)
  totals <- rowSums(pred$y)
  total_draws[[name]] <- totals
  scenario_tables[[name]] <- data.frame(
    estado = "Pendiente de validacion",
    escenario = name,
    viviendas = sum(pred$viviendas),
    segmentos = nrow(destino),
    censado_mismo_destino = sum(destino$n_per),
    t(summarize_total(totals)),
    check.names = FALSE
  )
  intervals <- t(apply(pred$y, 2, quantile, c(0.025, 0.5, 0.975)))
  intensity_intervals <- t(apply(pred$lambda, 2, quantile, c(0.025, 0.5, 0.975)))
  segment_tables[[name]] <- data.frame(
    escenario = name,
    destino[c("edid", "fid", "cut", "zona", "n_per", "n_vpo", "apto_modelo")],
    cobertura_asumida = pred$cobertura,
    V_tot = pred$viviendas,
    media = colMeans(pred$y),
    inferior = intervals[, 1],
    mediana = intervals[, 2],
    superior = intervals[, 3],
    intensidad_inferior = intensity_intervals[, 1],
    intensidad_mediana = intensity_intervals[, 2],
    intensidad_superior = intensity_intervals[, 3]
  )
  groups <- sort(unique(as.character(destino$cut)))
  commune_draws[[name]] <- sapply(groups, function(g)
    rowSums(pred$y[, as.character(destino$cut) == g, drop = FALSE]))
  commune_tables[[name]] <- do.call(rbind, lapply(groups, function(g) {
    columns <- as.character(destino$cut) == g
    data.frame(
      escenario = name,
      comuna = g,
      segmentos = sum(columns),
      censado_mismo_destino = sum(destino$n_per[columns]),
      t(summarize_total(commune_draws[[name]][, g])),
      check.names = FALSE
    )
  }))
  if (name == "base") {
    base_y <- pred$y
    base_lambda <- pred$lambda
  }
}

scenario_summary <- do.call(rbind, scenario_tables)
segment_summary <- do.call(rbind, segment_tables)
commune_summary <- do.call(rbind, commune_tables)

scenario_summary$cambio_mediana_pct <- 100 *
  (scenario_summary$mediana / scenario_summary$mediana[scenario_summary$escenario == "base"] - 1)
print(scenario_summary, digits = 6)

## ---- prediccion-agregar ----
# 12101 es la comuna del ejemplo. Cambiar el código por otro presente en destino.
# Comparar los percentiles del total con el resultado de sumar límites por segmento.
# Cada fila comparte beta, gamma y escalas: sumar dentro de esa misma fila.

selected_commune <- "12101"
columns <- as.character(destino$cut) == selected_commune
stopifnot(any(columns))
selected_total <- rowSums(base_y[, columns, drop = FALSE])
print(summarize_total(selected_total))

# Identidad de control: sumar comunas por simulacion reproduce el total del destino.
stopifnot(isTRUE(all.equal(as.numeric(rowSums(commune_draws$base)), as.numeric(total_draws$base))))

# Los percentiles del total dependen de cómo varían juntos los segmentos.
aggregation_check <- data.frame(
  limite_inferior_suma_segmentos = sum(segment_summary$inferior[segment_summary$escenario == "base"]),
  limite_inferior_total_correcto = unname(quantile(total_draws$base, 0.025)),
  limite_superior_suma_segmentos = sum(segment_summary$superior[segment_summary$escenario == "base"]),
  limite_superior_total_correcto = unname(quantile(total_draws$base, 0.975)))
print(aggregation_check)

## ---- prediccion-figuras ----
# Comparar medianas e intervalos entre escenarios y frente al censo del mismo destino.
# El gráfico por comuna usa escala logarítmica en el eje del conteo.

g_scenarios <- ggplot(scenario_summary, aes(x = mediana, y = escenario)) +
  geom_segment(aes(x = inferior, xend = superior, yend = escenario), linewidth = 1, colour = "#1f4e79") +
  geom_point(size = 3, colour = "#1f4e79") +
  geom_vline(xintercept = sum(destino$n_per), linetype = "dashed", colour = "#c55a11") +
  scale_x_continuous(labels = spanish_labels) +
  scale_y_discrete(labels = c(base = "Base: exposición final", uniforme_095 = "Hipótesis: c = 0,95",
    focalizada_090 = "Hipótesis focalizada: c = 0,90")) +
  labs(x = "Conteo predictivo del destino: mediana e intervalo del 95 %", y = NULL,
    subtitle = "Universo del modelo final de Magallanes; línea naranja: conteo censal comparable") +
  theme(plot.margin = margin(8, 35, 8, 8))

ggsave("resultados/escenarios.png", g_scenarios, width = 12, height = 4.5, dpi = 120, bg = "white")

base_communes <- commune_summary[commune_summary$escenario == "base", ]

g_communes <- ggplot(base_communes, aes(y = reorder(comuna, mediana))) +
  geom_segment(aes(x = inferior, xend = superior, yend = reorder(comuna, mediana)), colour = "#1f4e79") +
  geom_point(aes(x = mediana), colour = "#1f4e79", size = 2.4) +
  geom_point(aes(x = censado_mismo_destino), colour = "#c55a11", shape = 4, size = 3) +
  scale_x_log10(labels = spanish_labels) + labs(x = "Personas en el destino de cada comuna (escala logarítmica)",
    y = "Comuna", subtitle = "Azul: predicción base; cruz naranja: censado en el mismo universo")

ggsave("resultados/comunas.png", g_communes, width = 11, height = 5.5, dpi = 120, bg = "white")

## ---- prediccion-guardar ----
# Los CSV sirven para revisar y compartir tablas. Los RDS conservan las réplicas
# para repetir las sumas y el resumen guarda semillas, firmas y universo.
example_ids <- unique(c(which.min(destino$n_vpo), order(destino$n_vpo)[ceiling(nrow(destino)/2)],
  which.max(destino$n_vpo)))

examples <- segment_tables$base[example_ids, ]
segment_summary$estado <- "Pendiente de validacion del muestreo"
commune_summary$estado <- "Pendiente de validacion del muestreo"

write.csv(segment_summary, "resultados/prediccion_segmentos.csv", row.names = FALSE)
write.csv(commune_summary, "resultados/totales_comuna.csv", row.names = FALSE)
write.csv(scenario_summary, "resultados/escenarios.csv", row.names = FALSE)

saveRDS(total_draws, "resultados/totales_draws.rds")
saveRDS(commune_draws, "resultados/comunas_draws.rds")
saveRDS(list(escenarios = scenario_summary, comunas = commune_summary,
  ejemplos_segmentos = examples, agregacion = aggregation_check,
  destino_segmentos = nrow(destino), ajuste_segmentos = nrow(fit_data$observed),
  censado_destino = sum(destino$n_per), cap_pct = 100 * mean(density_raw > delta_max),
  estado = "Pendiente de validacion", configuracion = list(n_draws = length(ids), ids = ids, delta_max = delta_max,
    seed_density = 2621, seed_uniforms = 2622, stan_hash = config$stan_hash,
    columns_X = config$columns_X, groups = config$groups,
    source_chains = "../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds",
    chains_hash = unname(tools::md5sum("../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")),
    target_hash = unname(tools::md5sum("data/destino_magallanes.rds")))),
  "resultados/resumen.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")

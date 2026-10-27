# D4S2: prediccion desde las cadenas de D4S1; no ejecuta Stan.
# Ejecutar por bloques desde la carpeta Dia4_S2_Prediccion.

## ---- prediccion-paquetes ----
library(posterior)
library(ggplot2)
theme_set(theme_minimal(base_size = 14))
dir.create("resultados", showWarnings = FALSE)

## ---- prediccion-insumos ----
fit_data <- readRDS("../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")
diagnostics <- readRDS("../Dia4_S1_Diagnostico/resultados/resumen.rds")
target <- readRDS("data/destino_magallanes.rds")
destino <- target$destino
X_tot <- target$X_tot
Z_tot <- target$Z_tot
config <- fit_data$configuration
stopifnot(!diagnostics$sampler_review, identical(config, diagnostics$configuration),
  identical(target$stan_hash, config$stan_hash),
  identical(colnames(X_tot), config$columns_X),
  identical(colnames(Z_tot), config$groups),
  !anyDuplicated(destino[c("edid", "fid")]),
  nrow(destino) == nrow(X_tot), nrow(destino) == nrow(Z_tot),
  all(is.finite(X_tot)), all(is.finite(Z_tot)),
  all(rowSums(Z_tot) == 1), all(Z_tot %in% c(0,1)),
  all(as.character(destino$cut) %in% config$groups),
  identical(as.integer(max.col(Z_tot, ties.method = "first")),
    as.integer(match(as.character(destino$cut), config$groups))),
  all(destino$n_vpo > 0))
# Verificar que los aptos del ajuste esten representados sin cambiar sus datos.
key <- function(d) paste(d$edid, d$fid, sep = "\r")
idx <- match(key(fit_data$observed), key(destino))
stopifnot(!anyNA(idx),
  identical(as.numeric(destino$n_per[idx]), as.numeric(fit_data$observed$n_per)),
  identical(as.numeric(destino$n_vpo[idx]), as.numeric(fit_data$observed$n_vpo)),
  isTRUE(all.equal(unname(X_tot[idx, ]), unname(fit_data$X_obs))),
  isTRUE(all.equal(unname(Z_tot[idx, ]), unname(fit_data$Z_obs))))
cat("Destino: subconjunto sin colectivas y con VPO positiva. No es toda la region.\n")

## ---- prediccion-draws ----
# Seleccion espaciada para limitar memoria, solo para prediccion.
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
# Numeros aleatorios comunes: facilitan comparar exposiciones, no son efectos causales.
set.seed(2622)
uniforms <- matrix(runif(length(density)), nrow(density))

## ---- prediccion-funcion ----
predecir <- function(cobertura) {
  stopifnot(length(cobertura) == nrow(destino), all(is.finite(cobertura)),
            all(cobertura > 0 & cobertura <= 1))
  viviendas <- round(destino$n_vpo / cobertura)
  lambda <- pmin(pmax(sweep(density, 2, viviendas, "*"), 1), 1e7)
  # qpois transforma los mismos uniformes en conteos Poisson por escenario.
  y <- matrix(qpois(as.vector(uniforms), as.vector(lambda)), nrow(lambda))
  stopifnot(all(is.finite(y)), all(y >= 0))
  list(y = y, lambda = lambda, viviendas = viviendas, cobertura = cobertura)
}

## ---- prediccion-escenarios ----
# Son supuestos didacticos fijos, no tasas de cobertura estimadas.
c_base <- rep(1, nrow(destino))
c_uniforme <- rep(0.95, nrow(destino))
c_focalizada <- ifelse(as.character(destino$clasificacion_cobertura) == "Incompleta", 0.90, 1)
coverage_scenarios <- list(base = c_base, uniforme_095 = c_uniforme,
                           focalizada_090 = c_focalizada)

## ---- prediccion-llaves ----
# Ejemplo de archivo externo desordenado: alinear por edid/fid, nunca por posicion.
alinear_cobertura <- function(tabla) {
  stopifnot(is.data.frame(tabla), all(c("edid", "fid", "cobertura") %in% names(tabla)),
    is.numeric(tabla$cobertura), !anyDuplicated(tabla[c("edid", "fid")]),
    nrow(tabla) == nrow(destino))
  reference <- data.frame(destino[c("edid", "fid")], orden = seq_len(nrow(destino)))
  joined <- merge(reference, tabla, by = c("edid", "fid"), all.x = TRUE, sort = FALSE)
  joined <- joined[order(joined$orden), ]
  stopifnot(nrow(joined) == nrow(destino), !anyNA(joined$cobertura),
            all(is.finite(joined$cobertura)),
            all(joined$cobertura > 0 & joined$cobertura <= 1))
  joined$cobertura
}
set.seed(2623)
coverage_table <- data.frame(destino[c("edid", "fid")], cobertura = c_focalizada)
coverage_table <- coverage_table[sample(nrow(coverage_table)), ]
stopifnot(identical(as.numeric(alinear_cobertura(coverage_table)), as.numeric(c_focalizada)))

## ---- prediccion-resumenes ----
summarize_total <- function(x) {
  c(media = mean(x), mediana = median(x),
    inferior = unname(quantile(x, 0.025)), superior = unname(quantile(x, 0.975)))
}
scenario_tables <- segment_tables <- commune_tables <- list()
total_draws <- commune_draws <- list()
base_y <- base_lambda <- NULL
for (name in names(coverage_scenarios)) {
  pred <- predecir(coverage_scenarios[[name]])
  totals <- rowSums(pred$y)
  total_draws[[name]] <- totals
  scenario_tables[[name]] <- data.frame(escenario = name,
    viviendas = sum(pred$viviendas), segmentos = nrow(destino),
    censado_mismo_destino = sum(destino$n_per), t(summarize_total(totals)), check.names = FALSE)
  intervals <- t(apply(pred$y, 2, quantile, c(0.025, 0.5, 0.975)))
  intensity_intervals <- t(apply(pred$lambda, 2, quantile, c(0.025, 0.5, 0.975)))
  segment_tables[[name]] <- data.frame(escenario = name,
    destino[c("edid", "fid", "cut", "zona", "n_per", "n_vpo", "apto_modelo")],
    cobertura_asumida = pred$cobertura, V_tot = pred$viviendas,
    media = colMeans(pred$y), inferior = intervals[, 1], mediana = intervals[, 2],
    superior = intervals[, 3], intensidad_inferior = intensity_intervals[, 1],
    intensidad_mediana = intensity_intervals[, 2], intensidad_superior = intensity_intervals[, 3])
  groups <- sort(unique(as.character(destino$cut)))
  commune_draws[[name]] <- sapply(groups, function(g)
    rowSums(pred$y[, as.character(destino$cut) == g, drop = FALSE]))
  commune_tables[[name]] <- do.call(rbind, lapply(groups, function(g) {
    columns <- as.character(destino$cut) == g
    data.frame(escenario = name, comuna = g, segmentos = sum(columns),
      censado_mismo_destino = sum(destino$n_per[columns]),
      t(summarize_total(commune_draws[[name]][, g])), check.names = FALSE)
  }))
  if (name == "base") {base_y <- pred$y; base_lambda <- pred$lambda}
}
scenario_summary <- do.call(rbind, scenario_tables)
segment_summary <- do.call(rbind, segment_tables)
commune_summary <- do.call(rbind, commune_tables)
scenario_summary$cambio_mediana_pct <- 100 *
  (scenario_summary$mediana / scenario_summary$mediana[scenario_summary$escenario == "base"] - 1)
print(scenario_summary, digits = 6)

## ---- prediccion-agregar ----
# Cada fila comparte beta, gamma y escalas: sumar dentro de esa misma fila.
selected_commune <- "12101"
columns <- as.character(destino$cut) == selected_commune
stopifnot(any(columns))
selected_total <- rowSums(base_y[, columns, drop = FALSE])
print(summarize_total(selected_total))
# Identidad de control: sumar comunas por simulacion reproduce el total del destino.
stopifnot(isTRUE(all.equal(as.numeric(rowSums(commune_draws$base)), as.numeric(total_draws$base))))
# La media si es aditiva; medianas y limites no tienen por que serlo.
aggregation_check <- data.frame(
  limite_inferior_suma_segmentos = sum(segment_summary$inferior[segment_summary$escenario == "base"]),
  limite_inferior_total_correcto = unname(quantile(total_draws$base, 0.025)),
  limite_superior_suma_segmentos = sum(segment_summary$superior[segment_summary$escenario == "base"]),
  limite_superior_total_correcto = unname(quantile(total_draws$base, 0.975)))
print(aggregation_check)

## ---- prediccion-figuras ----
spanish_labels <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
g_scenarios <- ggplot(scenario_summary, aes(x = mediana, y = escenario)) +
  geom_segment(aes(x = inferior, xend = superior, yend = escenario), linewidth = 1, colour = "#1f4e79") +
  geom_point(size = 3, colour = "#1f4e79") +
  geom_vline(xintercept = sum(destino$n_per), linetype = "dashed", colour = "#c55a11") +
  scale_x_continuous(labels = spanish_labels) +
  scale_y_discrete(labels = c(base = "Base: c = 1", uniforme_095 = "Hipótesis: c = 0,95",
    focalizada_090 = "Hipótesis focalizada: c = 0,90")) +
  labs(x = "Conteo predictivo del destino: mediana e intervalo del 95 %", y = NULL,
    subtitle = "Mismo subconjunto de Magallanes; línea naranja: conteo censal comparable") +
  theme(plot.margin = margin(8, 35, 8, 8))
ggsave("resultados/escenarios.png", g_scenarios, width = 12, height = 4.5, dpi = 120, bg = "white")
base_communes <- commune_summary[commune_summary$escenario == "base", ]
g_communes <- ggplot(base_communes, aes(y = reorder(comuna, mediana))) +
  geom_segment(aes(x = inferior, xend = superior, yend = reorder(comuna, mediana)), colour = "#1f4e79") +
  geom_point(aes(x = mediana), colour = "#1f4e79", size = 2.4) +
  geom_point(aes(x = censado_mismo_destino), colour = "#c55a11", shape = 4, size = 3) +
  scale_x_log10(labels = spanish_labels) + labs(x = "Personas en el destino de cada comuna (escala logarítmica)",
    y = "Comuna", subtitle = "Azul: predicción base; cruz naranja: censado en el mismo subconjunto")
ggsave("resultados/comunas.png", g_communes, width = 11, height = 5.5, dpi = 120, bg = "white")

## ---- prediccion-guardar ----
example_ids <- unique(c(which.min(destino$n_vpo), order(destino$n_vpo)[ceiling(nrow(destino)/2)],
  which.max(destino$n_vpo)))
examples <- segment_tables$base[example_ids, ]
write.csv(segment_summary, "resultados/prediccion_segmentos.csv", row.names = FALSE)
write.csv(commune_summary, "resultados/totales_comuna.csv", row.names = FALSE)
write.csv(scenario_summary, "resultados/escenarios.csv", row.names = FALSE)
saveRDS(total_draws, "resultados/totales_draws.rds")
saveRDS(commune_draws, "resultados/comunas_draws.rds")
saveRDS(list(escenarios = scenario_summary, comunas = commune_summary,
  ejemplos_segmentos = examples, agregacion = aggregation_check,
  destino_segmentos = nrow(destino), ajuste_segmentos = nrow(fit_data$observed),
  censado_destino = sum(destino$n_per), cap_pct = 100 * mean(density_raw > delta_max),
  configuracion = list(n_draws = length(ids), ids = ids, delta_max = delta_max,
    seed_density = 2621, seed_uniforms = 2622, stan_hash = config$stan_hash,
    columns_X = config$columns_X, groups = config$groups,
    source_chains = "../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds",
    chains_hash = unname(tools::md5sum("../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")),
    target_hash = unname(tools::md5sum("data/destino_magallanes.rds")))),
  "resultados/resumen.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")

# D4S3: practica integradora. Ejecutar por bloques desde esta carpeta.
# Usa las cadenas de D4S1 y el destino de D4S2. No ajusta otro modelo Stan.

## ---- practica-paquetes ----
library(posterior)
library(bayesplot)
library(ggplot2)
theme_set(theme_minimal(base_size = 14))
color_scheme_set("viridis")
dir.create("resultados", showWarnings = FALSE)

## ---- practica-insumos ----
fit_data <- readRDS("../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")
target <- readRDS("../Dia4_S2_Prediccion/data/destino_magallanes.rds")
config <- fit_data$configuration
observado <- fit_data$observed
destino <- target$destino
X_tot <- target$X_tot
Z_tot <- target$Z_tot
stopifnot(identical(target$stan_hash, config$stan_hash),
  identical(target$source_configuration, config),
  identical(colnames(X_tot), config$columns_X), identical(colnames(Z_tot), config$groups),
  !anyDuplicated(destino[c("edid", "fid")]),
  all(is.finite(X_tot)), all(is.finite(Z_tot)),
  all(Z_tot %in% c(0,1)), all(rowSums(Z_tot) == 1),
  all(as.character(destino$cut) %in% config$group_labels),
  identical(as.integer(max.col(Z_tot, ties.method = "first")),
            as.integer(match(as.character(destino$cut), config$group_labels))),
  nrow(destino) == nrow(X_tot), nrow(destino) == nrow(Z_tot),
  all(is.finite(destino$n_vpo)), all(destino$n_vpo >= 0), all(is.finite(destino$V_tot)), all(destino$V_tot >= 0))
key <- function(d) paste(d$edid, d$fid, sep = "\r")
idx <- match(key(observado), key(destino))
stopifnot(!anyNA(idx), identical(as.numeric(observado$n_per), as.numeric(destino$n_per[idx])),
  identical(as.numeric(observado$n_vpo), as.numeric(destino$n_vpo[idx])),
  isTRUE(all.equal(unname(fit_data$X_obs), unname(X_tot[idx, ]))),
  isTRUE(all.equal(unname(fit_data$Z_obs), unname(Z_tot[idx, ]))))
print(data.frame(conjunto = c("Aptos del ajuste", "Destino de prediccion"),
  segmentos = c(nrow(observado), nrow(destino))))

## ---- practica-diagnostico ----
draws <- posterior::as_draws_array(fit_data$draws)
diagnostico <- posterior::summarise_draws(draws,
  "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")
diagnostico$revisar <- with(diagnostico,
  !is.finite(rhat) | rhat >= 1.01 | ess_bulk < 400 | ess_tail < 400)
nuts <- do.call(rbind, lapply(seq_along(fit_data$sampler), function(k) {
  s <- fit_data$sampler[[k]]
  energy <- s[, "energy__"]
  data.frame(cadena = k, divergencias = sum(s[, "divergent__"]),
    limite_profundidad = sum(s[, "treedepth__"] >= config$max_treedepth),
    ebfmi = sum(diff(energy)^2) / sum((energy - mean(energy))^2))
}))
nuts$revisar <- with(nuts, divergencias > 0 | limite_profundidad > 0 |
  !is.finite(ebfmi) | ebfmi < 0.3)
calculo_adecuado <- !any(diagnostico$revisar) && !any(nuts$revisar)
print(as.data.frame(diagnostico), digits = 4)
print(nuts, digits = 4)

## ---- practica-trazas ----
peores <- diagnostico$variable[order(diagnostico$rhat, decreasing = TRUE)][1:2]
g_trace <- mcmc_trace(draws, pars = peores, facet_args = list(ncol = 1)) +
  labs(x = "Iteración posterior", subtitle = "Revisar las variables de mayor R-hat junto a sus ESS") +
  guides(colour = guide_legend(title = "Cadena"), fill = guide_legend(title = "Cadena"))
ggsave("resultados/trazas.png", g_trace, width = 12, height = 5.5, dpi = 120, bg = "white")

## ---- practica-replicas ----
# Replica mixta de los aptos: nueva densidad de segmento, sin techo predictivo.
all_draws <- posterior::as_draws_matrix(draws)
beta_names <- paste0("beta[", seq_len(ncol(X_tot)), "]")
gamma_names <- paste0("gamma[", seq_len(ncol(Z_tot)), "]")
ids_check <- unique(round(seq(1, nrow(all_draws), length.out = min(500, nrow(all_draws)))))
beta_check <- all_draws[ids_check, beta_names, drop = FALSE]
gamma_check <- all_draws[ids_check, gamma_names, drop = FALSE]
sigma_check <- as.numeric(all_draws[ids_check, "sigma"])
eta_check <- beta_check %*% t(fit_data$X_obs) + gamma_check %*% t(fit_data$Z_obs)
set.seed(2612)
density_check <- exp(eta_check + matrix(rnorm(length(eta_check)), nrow(eta_check)) * sigma_check)
lambda_check <- pmin(pmax(sweep(density_check, 2, observado$n_vpo, "*"), 1), 1e7)
y_rep <- matrix(rpois(length(lambda_check), lambda_check), nrow(lambda_check))

## ---- practica-chequeos ----
# Cambie "zona" por "cut" y repita este bloque para comparar comunas.
variable_chequeo <- "zona"
stopifnot(variable_chequeo %in% c("zona", "cut"))
grupo <- as.character(observado[[variable_chequeo]])
chequeo <- do.call(rbind, lapply(sort(unique(grupo)), function(g) {
  columns <- grupo == g
  tasas <- rowSums(y_rep[, columns, drop = FALSE]) / sum(observado$n_vpo[columns])
  data.frame(grupo = g, segmentos = sum(columns),
    observado = sum(observado$n_per[columns]) / sum(observado$n_vpo[columns]),
    mediana_replica = median(tasas), inferior = unname(quantile(tasas, 0.025)),
    superior = unname(quantile(tasas, 0.975)))
}))
variacion_rep <- apply(sweep(y_rep, 2, observado$n_vpo, "/"), 1, sd)
variacion <- data.frame(observado = sd(observado$n_per / observado$n_vpo),
  inferior = unname(quantile(variacion_rep, 0.025)), mediana_replica = median(variacion_rep),
  superior = unname(quantile(variacion_rep, 0.975)))
print(chequeo)
print(variacion)

## ---- practica-posterior-predictiva ----
if (!calculo_adecuado) warning("DEMOSTRACION DE CALCULO: ajuste final sin convergencia validada; no usar para inferencia.")
ids_pred <- unique(round(seq(1, nrow(all_draws), length.out = min(1000, nrow(all_draws)))))
beta <- all_draws[ids_pred, beta_names, drop = FALSE]
gamma <- all_draws[ids_pred, gamma_names, drop = FALSE]
sigma <- as.numeric(all_draws[ids_pred, "sigma"])
eta <- beta %*% t(X_tot) + gamma %*% t(Z_tot)
set.seed(2621)
density <- pmin(exp(eta + matrix(rnorm(length(eta)), nrow(eta)) * sigma), config$delta_max)
set.seed(2622)
uniforms <- matrix(runif(length(density)), nrow(density))
predecir <- function(cobertura) {
  stopifnot(length(cobertura) == nrow(destino), all(is.finite(cobertura)),
    all(cobertura > 0 & cobertura <= 1))
  viviendas <- round(destino$V_tot / cobertura)
  lambda <- pmin(pmax(sweep(density, 2, viviendas, "*"), 1), 1e7)
  y <- matrix(qpois(as.vector(uniforms), as.vector(lambda)), nrow(lambda))
  stopifnot(all(is.finite(y)), all(y >= 0))
  list(y = y, viviendas = viviendas, cobertura = cobertura)
}

## ---- practica-eleccion ----
# Cambie la comuna y pruebe otro valor de cobertura. Son supuestos hipoteticos.
comuna_elegida <- "12101"
coberturas_probar <- c(0.90, 0.85)
stopifnot(length(comuna_elegida) == 1, comuna_elegida %in% as.character(destino$cut),
  length(coberturas_probar) >= 1, !anyDuplicated(coberturas_probar),
  all(is.finite(coberturas_probar)), all(coberturas_probar > 0 & coberturas_probar <= 1))
afectados <- as.character(destino$cut) == comuna_elegida &
  as.character(destino$clasificacion_cobertura) == "Incompleta"
if (!any(afectados)) stop("La comuna no tiene segmentos Incompleta en el destino. Revise la tabla y elija otra.")
escenarios_c <- list(base = rep(1, nrow(destino)))
for (c in coberturas_probar) {
  cobertura <- rep(1, nrow(destino))
  cobertura[afectados] <- c
  escenarios_c[[paste0("c_", format(c, trim = TRUE, scientific = FALSE))]] <- cobertura
}
print(data.frame(comuna = comuna_elegida, segmentos_del_grupo = sum(afectados)))

## ---- practica-resultados ----
resumir <- function(x) c(media = mean(x), mediana = median(x),
  inferior = unname(quantile(x, 0.025)), superior = unname(quantile(x, 0.975)))
columnas_comuna <- as.character(destino$cut) == comuna_elegida
tabla_destino <- tabla_comuna <- segmentos <- list()
totales_draws <- comunas_draws <- list()
for (nombre in names(escenarios_c)) {
  p <- predecir(escenarios_c[[nombre]])
  total <- rowSums(p$y)
  total_comuna <- rowSums(p$y[, columnas_comuna, drop = FALSE])
  totales_draws[[nombre]] <- total
  groups <- sort(unique(as.character(destino$cut)))
  comunas_draws[[nombre]] <- sapply(groups, function(g)
    rowSums(p$y[, as.character(destino$cut) == g, drop = FALSE]))
  stopifnot(identical(as.numeric(rowSums(comunas_draws[[nombre]])), as.numeric(total)))
  tabla_destino[[nombre]] <- data.frame(estado = "Pendiente de validacion", escenario = nombre,
    viviendas = sum(p$viviendas), segmentos_modificados = sum(p$cobertura != 1),
    censado_mismo_destino = sum(destino$n_per), t(resumir(total)), check.names = FALSE)
  tabla_comuna[[nombre]] <- data.frame(escenario = nombre, comuna = comuna_elegida,
    segmentos = sum(columnas_comuna), censado_mismo_destino = sum(destino$n_per[columnas_comuna]),
    t(resumir(total_comuna)), check.names = FALSE)
  intervals <- t(apply(p$y, 2, quantile, c(0.025, 0.5, 0.975)))
  segmentos[[nombre]] <- data.frame(escenario = nombre,
    destino[c("edid", "fid", "cut", "zona", "n_per", "n_vpo", "apto_modelo")],
    cobertura_asumida = p$cobertura, V_tot = p$viviendas, media = colMeans(p$y),
    inferior = intervals[, 1], mediana = intervals[, 2], superior = intervals[, 3])
}
tabla_destino <- do.call(rbind, tabla_destino)
tabla_comuna <- do.call(rbind, tabla_comuna)
segmentos <- do.call(rbind, segmentos)
tabla_comuna$cambio_mediana_pct <- 100 * (tabla_comuna$mediana / tabla_comuna$mediana[1] - 1)
print(tabla_destino)
print(tabla_comuna)

## ---- practica-agregacion ----
# Verificar una comuna diferente, sin cambiar el escenario base.
comuna_para_resumir <- "12401"
stopifnot(comuna_para_resumir %in% colnames(comunas_draws$base))
resumen_comuna <- resumir(comunas_draws$base[, comuna_para_resumir])
print(resumen_comuna)
# No sumar limites ni medianas de segmentos para fabricar un intervalo de comuna.

## ---- practica-figuras ----
g_check <- ggplot(chequeo, aes(y = grupo)) +
  geom_segment(aes(x = inferior, xend = superior, yend = grupo), colour = "#1f4e79", linewidth = 1) +
  geom_point(aes(x = mediana_replica), colour = "#1f4e79", size = 3) +
  geom_point(aes(x = observado), colour = "#c55a11", shape = 4, size = 3.5) +
  labs(x = "Personas por vivienda en aptos", y = variable_chequeo,
    subtitle = "Azul: réplicas mixtas e intervalo 95 %; cruz naranja: observado")
ggsave("resultados/chequeo.png", g_check, width = 11, height = 5, dpi = 120, bg = "white")
spanish_labels <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
g_scenarios <- ggplot(tabla_comuna, aes(y = escenario, x = mediana)) +
  geom_segment(aes(x = inferior, xend = superior, yend = escenario), colour = "#1f4e79", linewidth = 1) +
  geom_point(colour = "#1f4e79", size = 3) +
  geom_vline(xintercept = sum(destino$n_per[columnas_comuna]), colour = "#c55a11", linetype = "dashed") +
  scale_x_continuous(labels = spanish_labels) +
  labs(x = "Conteo predictivo: mediana e intervalo del 95 %", y = NULL,
    subtitle = paste("Comuna", comuna_elegida, "solo en el destino definido; naranja: censado comparable")) +
  theme(plot.margin = margin(8, 35, 8, 8))
ggsave("resultados/escenarios_comuna.png", g_scenarios, width = 12, height = 4.5, dpi = 120, bg = "white")

## ---- practica-guardar ----
tabla_comuna$estado <- "Pendiente de validacion del muestreo"
segmentos$estado <- "Pendiente de validacion del muestreo"
write.csv(as.data.frame(diagnostico), "resultados/diagnostico_parametros.csv", row.names = FALSE)
write.csv(nuts, "resultados/diagnostico_cadenas.csv", row.names = FALSE)
write.csv(chequeo, "resultados/chequeo_grupos.csv", row.names = FALSE)
write.csv(tabla_destino, "resultados/escenarios_destino.csv", row.names = FALSE)
write.csv(tabla_comuna, "resultados/escenarios_comuna.csv", row.names = FALSE)
write.csv(segmentos, "resultados/prediccion_segmentos.csv", row.names = FALSE)
saveRDS(totales_draws, "resultados/totales_draws.rds")
saveRDS(comunas_draws, "resultados/comunas_draws.rds")
saveRDS(list(diagnostico = as.data.frame(diagnostico), nuts = nuts,
  estado = "Pendiente de validacion", calculo_adecuado = calculo_adecuado, chequeo = chequeo, variacion = variacion,
  tabla_destino = tabla_destino, tabla_comuna = tabla_comuna,
  comuna_elegida = comuna_elegida, coberturas_probar = coberturas_probar,
  segmentos_grupo = sum(afectados), variable_chequeo = variable_chequeo,
  n_observado = nrow(observado), n_destino = nrow(destino),
  configuracion = list(ids_check = ids_check, ids_pred = ids_pred,
    seed_ppc = 2612, seed_density = 2621, seed_uniforms = 2622,
    delta_max = config$delta_max, columns_X = config$columns_X, groups = config$groups,
    chains_hash = unname(tools::md5sum("../Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")),
    target_hash = unname(tools::md5sum("../Dia4_S2_Prediccion/data/destino_magallanes.rds")))),
  "resultados/resumen.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")

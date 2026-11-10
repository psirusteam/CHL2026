# Ejecutar desde v02_2026-10-21. Comprueba datos y cadenas originales.
library(posterior)
source_dir <- "../../../output/04_seleccion_covariables_generalizado/cut_zona_fija/ajuste"
report <- list()
for (region in c("11", "12")) {
 x <- readRDS(file.path("Modelo_final", paste0("region_", region, ".rds")))
 folder <- file.path(source_dir, paste0("region_", region))
 original_data <- readRDS(file.path(folder, "datos_stan.rds"))
 original_chains <- readRDS(file.path(folder, paste0("cadenas_region_", region, ".rds")))
 stopifnot(identical(x$data_stan, original_data), length(x$covariates) == 40L,
   ncol(x$X_obs) == 43L, dim(x$draws)[1] == 500L, dim(x$draws)[2] == 4L,
   all(x$covariates %in% names(x$segments)), !anyDuplicated(x$segments[c("edid", "fid")]),
   identical(as.numeric(x$segments$V_tot), as.numeric(original_data$V_tot)),
   identical(x$X_obs, original_data$X_obs), identical(x$Z_obs, original_data$Z_obs),
   identical(as.numeric(x$observed$n_per), as.numeric(original_data$Y_obs)))
 mat <- as.matrix(posterior::as_draws_matrix(posterior::as_draws_array(x$draws)))
 vars <- colnames(mat)
 chain_order <- order(original_chains$.chain, original_chains$.iteration)
 stopifnot(isTRUE(all.equal(as.numeric(mat), as.numeric(as.matrix(original_chains[chain_order, vars])), check.attributes = FALSE)),
   isTRUE(all.equal(as.numeric(rowSums(x$Y_tot)), as.numeric(original_chains$Y_total_region[chain_order]))))
 hashes <- tools::md5sum(names(x$configuration$input_hashes))
 stopifnot(identical(unname(hashes), unname(x$configuration$input_hashes)))
 d <- posterior::summarise_draws(posterior::as_draws_array(x$draws), "rhat", "ess_bulk", "ess_tail")
 report[[region]] <- data.frame(region = region, covariables = length(x$covariates),
   P = original_data$P, Q = original_data$Q, D_obs = original_data$D_obs, D_tot = original_data$D_tot,
   rank_X = qr(original_data$X_obs)$rank, Rhat_max = max(d$rhat), ESS_bulk_min = min(d$ess_bulk),
   estado = "Pendiente de validacion del muestreo")
}
fit <- readRDS("Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")
target <- readRDS("Dia4_S2_Prediccion/data/destino_magallanes.rds")
x <- readRDS("Modelo_final/region_12.rds")
stopifnot(identical(fit$draws, x$draws), identical(fit$sampler, x$sampler),
  identical(target$X_tot, x$data_stan$X_tot), identical(target$Z_tot, x$data_stan$Z_tot),
  identical(target$V_tot, x$data_stan$V_tot), identical(target$destino, x$segments),
  identical(target$source_configuration, fit$configuration))
d <- readRDS("Dia4_S1_Diagnostico/resultados/resumen.rds")
p <- readRDS("Dia4_S2_Prediccion/resultados/resumen.rds")
q <- readRDS("Dia4_S3_Practica/resultados/resumen.rds")
stopifnot(d$sampler_review, !q$calculo_adecuado, p$destino_segmentos == nrow(target$destino), q$n_destino == nrow(target$destino))
a <- readRDS("Dia4_S2_Prediccion/resultados/totales_draws.rds")
b <- readRDS("Dia4_S3_Practica/resultados/totales_draws.rds")
stopifnot(identical(a$base, b$base))
report <- do.call(rbind, report)
write.csv(report, "Modelo_final/validacion_equivalencia.csv", row.names = FALSE)
print(report)
cat("Equivalencia y continuidad verificadas.\n")


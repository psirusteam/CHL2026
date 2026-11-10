# Ejecutar desde la carpeta v02_2026-10-21. Exporta, sin reestimar.
library(rstan)
library(posterior)
root <- normalizePath("../../..", winslash = "/")
source_dir <- file.path(root, "output/04_seleccion_covariables_generalizado/cut_zona_fija/ajuste")
stan_path <- file.path(root, "stan/modelo_poblacion_v02.stan")
stan_hash <- unname(tools::md5sum(stan_path))
covariates <- read.csv(file.path(source_dir, "covariables_modelo.csv"))$variable
stopifnot(length(covariates) == 40L)
file.copy(file.path(source_dir, "covariables_modelo.csv"), "Modelo_final/covariables_modelo.csv", overwrite = TRUE)
for (region in c("11", "12")) {
 folder <- file.path(source_dir, paste0("region_", region))
 session <- if (region == "11") "Dia3_S3_Modelo_Chile" else "Dia3_S4_Practica"
 label <- if (region == "11") "aysen" else "magallanes"
 data <- readRDS(file.path(folder, "datos_stan.rds"))
 identifiers <- readRDS(file.path(folder, "identificadores.rds"))
 groups <- readRDS(file.path(folder, "configuracion_grupos.rds"))
 configuration <- readRDS(file.path(folder, "configuracion_ajuste.rds"))
 stopifnot(identical(groups$covariables, covariates), stan_hash %in% groups$firma_insumos, data$P == 43L, groups$variable == "cut")
 segments <- readRDS(file.path(session, "data", paste0("segmentos_", label, ".rds")))
 key <- function(x) paste(x$edid, x$fid, sep = "\r")
 row_id <- match(key(identifiers), key(segments))
 stopifnot(!anyNA(row_id), !anyDuplicated(row_id), length(row_id) == nrow(segments))
 segments <- segments[row_id, ]
 for (name in covariates) segments[[name]] <- data$X_tot[, name]
 segments$V_tot <- data$V_tot
 stopifnot(identical(as.integer(segments$n_per[segments$apto_modelo]), data$Y_obs),
  identical(as.integer(segments$n_vpo[segments$apto_modelo]), data$V_obs),
  isTRUE(all.equal(unname(as.matrix(segments[segments$apto_modelo, covariates])), unname(data$X_obs[, covariates]))))
 saveRDS(segments, file.path(session, "data", paste0("segmentos_", label, ".rds")))
 cat("Leyendo ajuste implementado, region", region, "\n"); flush.console()
 fit <- readRDS(file.path(folder, paste0("ajuste_stan_region_", region, ".rds")))
 chains <- as.array(fit, pars = c("beta", "gamma", "sigma", "sigmaU"))
 sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
 y_tot <- as.matrix(fit, pars = "Y_tot")
 config <- list(region = label, seed_Stan = configuration$semilla, chains = configuration$cadenas,
  iter = configuration$iter, warmup = configuration$warmup, max_treedepth = configuration$max_treedepth,
  delta_max = data$delta_max, columns_X = colnames(data$X_obs), groups = colnames(data$Z_obs),
  group_labels = groups$niveles, covariates = covariates, stan_hash = stan_hash,
  source = file.path("output/04_seleccion_covariables_generalizado/cut_zona_fija/ajuste", paste0("region_", region)),
  input_hashes = tools::md5sum(file.path(folder, c("datos_stan.rds", "configuracion_grupos.rds", "configuracion_ajuste.rds"))),
  session = capture.output(sessionInfo()))
 saved <- list(draws = chains, sampler = sampler, observed = segments[segments$apto_modelo, ],
  X_obs = data$X_obs, Z_obs = data$Z_obs, configuration = config,
  data_stan = data, segments = segments, covariates = covariates,
  original_configuration = configuration, original_groups = groups, Y_tot = y_tot)
 saveRDS(saved, file.path("Modelo_final", paste0("region_", region, ".rds")))
 if (region == "12") {
  saved$Y_tot <- NULL
  saveRDS(saved, "Dia4_S1_Diagnostico/data/ajuste_magallanes.rds")
  saveRDS(list(destino = segments, X_tot = data$X_tot, Z_tot = data$Z_tot,
   V_tot = data$V_tot, stan_hash = stan_hash, source_configuration = config), "Dia4_S2_Prediccion/data/destino_magallanes.rds")
 }
 file.copy(stan_path, file.path(session, "stan/modelo_poblacion_v02.stan"), overwrite = TRUE)
 cat("Exportado", region, "P=",data$P,"Q=",data$Q,"D_tot=",data$D_tot,"\n")
 rm(fit, y_tot, saved); gc()
}

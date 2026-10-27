# Preparacion docente: regenerar cadenas completas para diagnostico.
# Ejecutar desde Dia4_S1_Diagnostico. No cambia salidas de D3S4.
# Se usa la misma especificacion, datos y configuracion de muestreo.
# Se fija ademas la semilla R de los valores iniciales.
library(rstan)
practice_dir <- normalizePath("../Dia3_S4_Practica", winslash = "/")
output_file <- normalizePath("data", winslash = "/")
output_file <- file.path(output_file, "ajuste_magallanes.rds")
lines <- readLines(file.path(practice_dir, "Rcode/Practica_Dia3.R"), encoding = "UTF-8")
end_line <- grep("^## ---- practica-resultados ----", lines)
stopifnot(length(end_line) == 1)
original_dir <- getwd()
setwd(practice_dir)
set.seed(2543)
source(textConnection(lines[seq_len(end_line - 1)]), encoding = "UTF-8")
# Preservar orden de iteraciones, cadenas y variables: sin mezclar ni adelgazar.
variables <- c("beta", "gamma", "sigma", "sigmaU")
saved <- list(
  draws = as.array(ajuste, pars = variables),
  sampler = rstan::get_sampler_params(ajuste, inc_warmup = FALSE),
  observed = aptos[c("edid", "fid", "cut", "zona", "n_per", "n_vpo")],
  X_obs = X_obs, Z_obs = Z_obs,
  configuration = list(region = "Magallanes", seed_R = 2543, seed_Stan = 2543,
    chains = 4, iter = 3000, warmup = 1500, adapt_delta = 0.99,
    max_treedepth = 12, delta_max = delta_max,
    columns_X = colnames(X_obs), groups = niveles_comuna,
    source_practice = "Dia3_S4_Practica/Rcode/Practica_Dia3.R",
    stan_hash = unname(tools::md5sum("stan/modelo_poblacion_v02.stan")),
    session = capture.output(sessionInfo())))
saveRDS(saved, output_file)
setwd(original_dir)
cat("Cadenas y diagnosticos guardados sin modificar los resultados de D3S4.\n")

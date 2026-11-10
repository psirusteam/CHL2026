# Utilidades compartidas de D3S3 a D4S3. Misma especificacion final.
#' Load a regional reference fit and check its model specification
#'
#' Read the compact package exported from the implemented Chile model.
#' Check that the local Stan file and selected design metadata match the
#' stored configuration. This function does not compile Stan or refit a model.
#'
#' @param region Regional code used in the package filename, for example
#'   "11" for Aysen or "12" for Magallanes. Use two digits to match the files.
#'
#' @return The complete list read from the RDS file, returned visibly.
#'   Components used in subsequent sessions include:
#'   * `draws`: post-warmup array, iteration x chain x parameter.
#'   * `sampler`: one matrix of NUTS diagnostics per chain.
#'   * `data_stan`: original counts, exposures and design matrices.
#'   * `segments` and `observed`: prediction and fitting records.
#'   * `covariates`: ordered names of the 40 standardized covariates.
#'   * `configuration`: Stan hash, design columns and sampling settings.
#'   * `Y_tot`: original predictive counts, draw x destination segment.
#'
#' @details
#' Paths are relative to the current working directory, which must be a
#' session folder such as Dia3_S3_Modelo_Chile. The package is read from
#' ../Modelo_final/region_<region>.rds; the local Stan file is read from stan/.
#'
#' Checks cover the exact Stan MD5 hash, 40 covariates, P = 43, the names and
#' order of X_obs columns, and the names and order of Z_tot columns.
#' A failed check raises an error through stopifnot(); an unreadable file
#' raises the corresponding readRDS() error. The list is never repaired.
#'
#' These checks do not establish convergence, validate population estimates,
#' or compare every matrix value or segment key. Use check_final_model() for
#' sampling diagnostics and Modelo_final/Validar_insumos.R, from the revision
#' root, for the broader comparison with the source project.
#'
#' @examples
#' # Run from Dia3_S3_Modelo_Chile or Dia3_S4_Practica:
#' # source("../Modelo_final/utilidades_modelo_final.R")
#' # final <- load_final_model("12")
#' # final$configuration$columns_X
#' 
load_final_model <- function(region) {
 x <- readRDS(file.path("Modelo_final", paste0("region_", region, ".rds")))
 stopifnot(identical(unname(tools::md5sum("stan/modelo_poblacion_v02.stan")), x$configuration$stan_hash),
  length(x$covariates) == 40L, x$data_stan$P == 43L,
  identical(colnames(x$data_stan$X_obs), x$configuration$columns_X),
  identical(colnames(x$data_stan$Z_tot), x$configuration$groups))
 x
}

#' Summarize parameter diagnostics and flag sampling issues
#'
#' Summarize the stored post-warmup draws and combine selected NUTS alerts
#' across chains. No iterations are discarded, reordered or thinned here.
#'
#' @param x A regional reference list, normally returned by load_final_model().
#'   Required components are `draws`, `sampler` and
#'   `configuration$max_treedepth`. Each sampler matrix must contain
#'   `divergent__` and `treedepth__` columns for post-warmup transitions.
#'
#' @return A visible list with two components:
#'   * `diagnostics`: posterior summary table, one row per stored parameter,
#'     with variable, mean, sd, rhat, ess_bulk and ess_tail.
#'   * `controls`: one-row data frame containing rhat_max, ess_bulk_min,
#'     ess_tail_min, divergencias, limite_profundidad and revisar.
#'
#' @details
#' Requires the posterior package, called through its namespace.
#' R-hat compares chain exploration; bulk ESS describes effective sample
#' size for the distribution body, while tail ESS describes tail precision.
#' Divergences and transitions reaching or exceeding max_treedepth are
#' counted over all supplied sampler matrices, not averaged across chains.
#'
#' The review flag is TRUE when the maximum R-hat is non-finite or >= 1.01,
#' either minimum ESS is < 400, or at least one divergence or maximum-depth
#' transition is recorded. Reaching maximum tree depth is an efficiency
#' alert; it does not by itself establish biased inference.
#'
#' The checks concern only stored parameters. The workshop reference contains
#' beta, gamma, sigma and sigmaU; it excludes individual latent densities.
#' This function does not compute E-BFMI, Monte Carlo standard errors or
#' predictive checks. A FALSE review flag is not evidence of model adequacy.
#'
#' Inputs are assumed to follow the reference package contract. Missing
#' columns cause errors; NA values are not removed and may propagate into
#' summaries or the review flag. No inputs or files are modified, and a TRUE
#' flag is returned rather than stopping downstream execution.
#'
#' @examples
#' # final <- load_final_model("12")
#' # checks <- check_final_model(final)
#' # checks$controls
#' # checks$diagnostics
check_final_model <- function(x) {
 d <- posterior::summarise_draws(posterior::as_draws_array(x$draws),
   "mean", "sd", "rhat", "ess_bulk", "ess_tail")
 s <- x$sampler
 controls <- data.frame(rhat_max = max(d$rhat), ess_bulk_min = min(d$ess_bulk),
  ess_tail_min = min(d$ess_tail),
  divergencias = sum(vapply(s, function(z) sum(z[, "divergent__"]), numeric(1))),
  limite_profundidad = sum(vapply(s, function(z) sum(z[, "treedepth__"] >= x$configuration$max_treedepth), numeric(1))))
 controls$revisar <- with(controls, !is.finite(rhat_max) | rhat_max >= 1.01 |
  ess_bulk_min < 400 | ess_tail_min < 400 | divergencias > 0 | limite_profundidad > 0)
 list(diagnostics = d, controls = controls)
}

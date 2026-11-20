# Ejecutar desde Dia4_S1_Diagnostico. Lee un ajuste regional a la vez.
# CHL2026_AJUSTES permite indicar otra ubicación de las exportaciones.
library(posterior)
library(ggplot2)

ruta_ajustes <- Sys.getenv("CHL2026_AJUSTES", unset = file.path(
  "../../../../PopulationModels/CHL2026/output",
  "04_seleccion_covariables_generalizado/cut_zona_fija/ajuste"
))
regiones <- sprintf("%02d", 1:16)
nombres <- c("Tarapacá", "Antofagasta", "Atacama", "Coquimbo", "Valparaíso",
  "O’Higgins", "Maule", "Biobío", "La Araucanía", "Los Lagos", "Aysén",
  "Magallanes", "Metropolitana", "Los Ríos", "Arica y Parinacota", "Ñuble")

# Los paquetes docentes de estas dos regiones incluyen los controles NUTS.
ruta_aysen <- Sys.getenv("CHL2026_AJUSTE_AYSEN",
  unset = "../.Version/19112026/Modelo_final/region_11.rds")
rutas_nuts <- c("11" = ruta_aysen, "12" = "data/ajuste_magallanes.rds")

chequear_ajuste <- function(matriz, datos, coef_fijos, coef_comuna, seed) {
  n <- datos$D_obs
  stopifnot(n > 1L, nrow(datos$X_obs) == n, nrow(datos$Z_obs) == n,
            length(datos$Y_obs) == n, length(datos$V_obs) == n,
            all(is.finite(datos$Y_obs)), all(datos$Y_obs >= 0),
            all(is.finite(datos$V_obs)), all(datos$V_obs > 0))
  grupos <- cbind(aldea = datos$X_obs[, "zona_aldea"],
                  urbana = datos$X_obs[, "zona_urbana"])
  stopifnot(all(grupos %in% c(0, 1)), all(rowSums(grupos) <= 1))
  zona <- ifelse(grupos[, "aldea"] == 1, "aldea",
                 ifelse(grupos[, "urbana"] == 1, "urbana", "rural"))
  niveles <- sort(unique(zona))
  ids <- unique(round(seq(1, nrow(matriz), length.out = min(500L, nrow(matriz)))))
  m <- length(ids)
  replicas <- matriz[ids, , drop = FALSE]
  beta <- replicas[, coef_fijos, drop = FALSE]
  comuna <- replicas[, coef_comuna, drop = FALSE]
  sigma <- replicas[, "sigma"]
  stopifnot(all(sigma > 0))
  conteos_zona <- matrix(0, m, length(niveles), dimnames = list(NULL, niveles))
  suma <- cuadrados <- numeric(m)
  seleccion_curvas <- unique(round(seq(1, m, length.out = min(50L, m))))
  replicas_curvas <- matrix(NA_real_, length(seleccion_curvas), n)

  # Acumular por bloques para evitar una matriz de 500 por todos los segmentos.
  set.seed(seed)
  for (inicio in seq.int(1L, n, by = 1024L)) {
    j <- inicio:min(inicio + 1023L, n)
    eta <- beta %*% t(datos$X_obs[j, , drop = FALSE]) +
      comuna %*% t(datos$Z_obs[j, , drop = FALSE])
    ruido <- matrix(rnorm(length(eta)), nrow = m) * sigma
    lambda <- exp(sweep(eta + ruido, 2, log(datos$V_obs[j]), "+"))
    lambda <- pmin(pmax(lambda, 1), 1e7)
    y_rep <- matrix(rpois(length(lambda), lambda), nrow = m)
    dens <- sweep(y_rep, 2, datos$V_obs[j], "/")
    replicas_curvas[, j] <- dens[seleccion_curvas, , drop = FALSE]
    suma <- suma + rowSums(dens)
    cuadrados <- cuadrados + rowSums(dens^2)
    for (z in niveles) {
      columnas <- zona[j] == z
      conteos_zona[, z] <- conteos_zona[, z] + rowSums(y_rep[, columnas, drop = FALSE])
    }
  }
  intervalo <- function(x) unname(quantile(x, c(0.025, 0.5, 0.975)))
  comparar <- function(medida, valor, simulado) {
    q <- intervalo(simulado)
    data.frame(medida = medida, observado = valor, inferior = q[1],
      mediana_replica = q[2], superior = q[3],
      dentro = valor >= q[1] && valor <= q[3])
  }
  por_zona <- lapply(niveles, function(z) {
    i <- zona == z
    viviendas <- sum(datos$V_obs[i])
    comparar(z, sum(datos$Y_obs[i]) / viviendas, conteos_zona[, z] / viviendas)
  })
  dispersion <- sqrt(pmax(0, (cuadrados - suma^2 / n) / (n - 1)))
  chequeos <- rbind(do.call(rbind, por_zona),
    comparar("dispersión", sd(datos$Y_obs / datos$V_obs), dispersion),
    comparar("densidad regional", sum(datos$Y_obs) / sum(datos$V_obs),
             rowSums(conteos_zona) / sum(datos$V_obs)))
  list(chequeos = chequeos, observados = datos$Y_obs / datos$V_obs,
       replicas = replicas_curvas)
}

diagnosticar_region <- function(codigo, nombre) {
  carpeta <- file.path(ruta_ajustes, paste0("region_", codigo))
  x <- readRDS(file.path(carpeta, paste0("cadenas_region_", codigo, ".rds")))
  cfg <- readRDS(file.path(carpeta, "configuracion_ajuste.rds"))
  datos <- readRDS(file.path(carpeta, "datos_stan.rds"))
  fijos <- paste0("beta[", seq_len(datos$P), "]")
  comunas <- paste0("gamma[", seq_len(datos$Q), "]")
  variables <- c(fijos, comunas, "sigma", "sigmaU")
  n_iter <- cfg$iter - cfg$warmup
  stopifnot(all(variables %in% names(x)), nrow(x) == n_iter * cfg$cadenas,
    identical(sort(unique(x$.chain)), seq_len(cfg$cadenas)),
    !anyDuplicated(x[c(".chain", ".iteration")]),
    all(x$.iteration %in% seq_len(n_iter)))
  x <- x[order(x$.chain, x$.iteration), ]
  matriz <- as.matrix(x[, variables])
  stopifnot(all(is.finite(matriz)))
  cadenas <- array(matriz, dim = c(n_iter, cfg$cadenas, length(variables)),
                  dimnames = list(NULL, NULL, variables))
  tabla <- as.data.frame(summarise_draws(as_draws_array(cadenas),
    "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  tabla$revisar <- with(tabla, !is.finite(rhat) | !is.finite(ess_bulk) |
    !is.finite(ess_tail) | rhat >= 1.01 | ess_bulk < 400 | ess_tail < 400)
  tabla$codigo <- codigo
  tabla$region <- nombre
  nuts <- NULL
  if (codigo %in% names(rutas_nuts) && file.exists(rutas_nuts[[codigo]])) {
    paquete <- readRDS(rutas_nuts[[codigo]])
    stopifnot(isTRUE(all.equal(as.vector(paquete$draws[, , variables]),
                               as.vector(cadenas), tolerance = 0)))
    nuts <- do.call(rbind, lapply(seq_along(paquete$sampler), function(k) {
      s <- paquete$sampler[[k]]
      energia <- s[, "energy__"]
      data.frame(codigo = codigo, region = nombre, cadena = k,
        divergencias = sum(s[, "divergent__"]),
        limite_profundidad = sum(s[, "treedepth__"] >= cfg$max_treedepth),
        ebfmi = sum(diff(energia)^2) / sum((energia - mean(energia))^2))
    }))
  }
  ajuste_ppc <- chequear_ajuste(matriz, datos, fijos, comunas, 2612L + as.integer(codigo))
  ppc <- ajuste_ppc$chequeos
  curva <- bayesplot::ppc_dens_overlay(ajuste_ppc$observados, ajuste_ppc$replicas) +
    coord_cartesian(xlim = c(0, 6)) +
    labs(title = nombre, x = "Personas por vivienda", y = NULL) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  ppc$codigo <- codigo
  ppc$region <- nombre
  peor <- function(x, fun) if (all(is.finite(x))) fun(x) else NA_real_
  resumen <- data.frame(codigo = codigo, region = nombre, aptos = datos$D_obs,
    cadenas = cfg$cadenas, iteraciones = n_iter,
    rhat_max = peor(tabla$rhat, max), ess_bulk_min = peor(tabla$ess_bulk, min),
    ess_tail_min = peor(tabla$ess_tail, min),
    parametros_con_alerta = sum(tabla$revisar), parametros = nrow(tabla),
    nuts_disponible = !is.null(nuts))
  list(resumen = resumen, parametros = tabla, nuts = nuts, ppc = ppc, curva = curva)
}

resultados <- lapply(seq_along(regiones), function(i) {
  cat("Revisando", regiones[i], nombres[i], "\n")
  flush.console()
  diagnosticar_region(regiones[i], nombres[i])
})
nacional <- lapply(c("resumen", "parametros", "nuts", "ppc"), function(campo) {
  do.call(rbind, lapply(resultados, function(x) x[[campo]]))
})
names(nacional) <- c("resumen", "parametros", "nuts", "ppc")
stopifnot(identical(as.character(nacional$resumen$codigo), regiones),
          length(unique(nacional$ppc$codigo)) == length(regiones))
dir.create("resultados", showWarnings = FALSE)
for (campo in names(nacional)) {
  write.csv(nacional[[campo]], file.path("resultados", paste0("nacional_", campo, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}
saveRDS(nacional, "resultados/diagnosticos_regionales.rds")

# Las dos mitades mantienen el mismo orden de regiones y las mismas escalas.
library(ggplot2)
tema_regional <- theme_minimal(base_size = 17) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))
ppc <- nacional$ppc[nacional$ppc$medida != "densidad regional", ]
ppc$medida <- factor(ppc$medida,
  levels = c("rural", "aldea", "urbana", "dispersión"),
  labels = c("Rural", "Aldea", "Urbana", "Dispersión"))
ppc$region <- factor(ppc$region, levels = rev(nombres))
for (grupo in list(1:8, 9:16)) {
  sufijo <- paste(range(regiones[grupo]), collapse = "-")
  # Las capas invisibles fijan un rango común entre los dos gráficos por medida.
  p <- ggplot(ppc[ppc$codigo %in% regiones[grupo], ], aes(y = region)) +
    geom_blank(data = ppc, aes(x = inferior), inherit.aes = FALSE) +
    geom_blank(data = ppc, aes(x = superior), inherit.aes = FALSE) +
    geom_segment(aes(x = inferior, xend = superior, yend = region),
                 colour = "#286394", linewidth = 1.1) +
    geom_point(aes(x = mediana_replica), colour = "#286394", size = 2.4) +
    geom_point(aes(x = observado), colour = "#a3410a", shape = 4, size = 3, stroke = 1.2) +
    facet_wrap(~ medida, ncol = 4, scales = "free_x") +
    labs(x = "Personas por vivienda (tasas por zona y dispersión entre segmentos)",
         y = NULL) + tema_regional
  ggsave(file.path("resultados", paste0("ppc_regiones_", sufijo, ".png")), p,
         width = 13, height = 4.5, dpi = 120, bg = "white")

  densidad <- nacional$ppc[nacional$ppc$medida == "densidad regional", ]
  densidad$region <- factor(densidad$region, levels = rev(nombres))
  p <- ggplot(densidad[densidad$codigo %in% regiones[grupo], ], aes(y = region)) +
    geom_segment(aes(x = inferior, xend = superior, yend = region),
                 colour = "#286394", linewidth = 1.1) +
    geom_point(aes(x = mediana_replica), colour = "#286394", size = 2.7) +
    geom_point(aes(x = observado), colour = "#a3410a", shape = 4,
               size = 3.4, stroke = 1.2) +
    coord_cartesian(xlim = range(densidad[c("inferior", "superior", "observado")])) +
    labs(x = "Personas por vivienda", y = NULL) + tema_regional
  ggsave(file.path("resultados", paste0("densidades_regiones_", sufijo, ".png")), p,
         width = 12, height = 4.5, dpi = 120, bg = "white")
}

# Guardar una figura por región y cuatro grupos de cuatro para la presentación.
dir.create("resultados/ppc_curvas", showWarnings = FALSE)
for (i in seq_along(regiones)) {
  ggsave(file.path("resultados/ppc_curvas", paste0("region_", regiones[i], ".png")),
         resultados[[i]]$curva, width = 7, height = 4, dpi = 150, bg = "white")
}
for (inicio in seq(1, 16, by = 4)) {
  grupo <- inicio:(inicio + 3L)
  curvas <- lapply(resultados[grupo], function(x) x$curva)
  figura <- patchwork::wrap_plots(curvas, ncol = 2)
  sufijo <- paste(range(regiones[grupo]), collapse = "-")
  ggsave(file.path("resultados", paste0("ppc_curvas_", sufijo, ".png")), figura,
         width = 12, height = 5.3, dpi = 135, bg = "white")
}

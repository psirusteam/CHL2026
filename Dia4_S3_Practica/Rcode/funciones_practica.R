# Funciones del ejercicio D4S3. Cargar desde el guion de la práctica.

#' Leer las cadenas y el destino de una región
#'
#' @param codigo Código regional de dos dígitos, entre 01 y 16.
#' @param directorio Carpeta local con un archivo region_XX.rds por región.
#' @return Lista con exposición, matrices de diseño, parámetros y totales de Stan.
leer_insumo_regional <- function(codigo, directorio = "data/regiones") {
  stopifnot(length(codigo) == 1L, codigo %in% sprintf("%02d", 1:16))
  ruta <- file.path(directorio, paste0("region_", codigo, ".rds"))
  if (!file.exists(ruta)) stop("Falta el insumo regional: ", ruta)
  x <- readRDS(ruta)
  n <- length(x$V)
  s <- length(x$sigma)
  stopifnot(identical(x$codigo, codigo), n > 0L, nrow(x$X) == n,
    length(x$grupo_comuna) == n, ncol(x$X) == ncol(x$beta),
    nrow(x$beta) == s, nrow(x$gamma) == s, length(x$total_stan) == s,
    all(x$grupo_comuna %in% seq_len(ncol(x$gamma))),
    all(is.finite(x$V)), all(x$V >= 0), all(is.finite(x$sigma)), all(x$sigma > 0),
    all(is.finite(x$total_stan)), all(x$total_stan >= 0),
    length(x$delta_max) == 1L, is.finite(x$delta_max), x$delta_max > 0,
    length(x$cadena) == s, length(x$iteracion) == s,
    !anyNA(x$cadena), !anyNA(x$iteracion),
    !anyDuplicated(data.frame(cadena = x$cadena, iteracion = x$iteracion)))
  x
}

#' Construir una combinación conjunta de réplicas independientes
#'
#' @param n Cantidad de réplicas por región.
#' @param codigos Vector de códigos regionales.
#' @param semilla Semilla para repetir el emparejamiento.
#' @return Matriz con una permutación por región y una fila por réplica conjunta.
#' @details Fija la semilla de R. Se reutiliza la matriz en todos los escenarios.
emparejar_regiones <- function(n, codigos, semilla) {
  set.seed(semilla)
  vapply(codigos, function(codigo) sample.int(n), integer(n))
}

#' Resumir un conteo predictivo o una diferencia por réplica
#'
#' @param x Vector con una suma o diferencia para cada réplica.
#' @return Vector con media, percentil 2,5, mediana y percentil 97,5.
#' @examples
#' resumir_replicas(c(100, 110, 120, 130, 140))
resumir_replicas <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.5, 0.975), type = 7, names = FALSE)
  c(media = mean(x), inferior = q[1L], mediana = q[2L], superior = q[3L])
}

#' Simular conteos de una región bajo tres coberturas
#'
#' @param x Paquete de una región leído con leer_insumo_regional().
#' @param n_replicas Cantidad de filas a tomar de las cadenas guardadas.
#' @param cobertura Vector nombrado base, uniforme y focalizada, en (0, 1].
#' @param zona_focal Zona donde se aplica la cobertura focalizada.
#' @param region_focal NA para aplicar el escenario en todas las regiones, o un código.
#' @param semilla Semilla regional de la simulación.
#' @return Lista con totales por réplica, totales por zona y tamaños del destino.
#' @details Se comparten las normales y los uniformes entre escenarios. Se fija
#'   la semilla de R y se procesan los segmentos en bloques para limitar memoria.
simular_coberturas <- function(x, n_replicas, cobertura, zona_focal, region_focal, semilla) {
  zonas <- c("Rural", "Aldea", "Urbana")
  stopifnot(length(n_replicas) == 1L, is.finite(n_replicas),
    n_replicas == floor(n_replicas), n_replicas >= 2L, n_replicas <= length(x$sigma),
    identical(names(cobertura), c("base", "uniforme", "focalizada")),
    all(is.finite(cobertura)), all(cobertura > 0 & cobertura <= 1), cobertura[["base"]] == 1,
    length(zona_focal) == 1L, zona_focal %in% zonas,
    length(region_focal) == 1L, is.na(region_focal) || region_focal %in% sprintf("%02d", 1:16))
  orden <- order(x$cadena, x$iteracion)
  filas <- orden[round(seq(1, length(orden), length.out = n_replicas))]
  beta <- x$beta[filas, , drop = FALSE]
  gamma <- x$gamma[filas, , drop = FALSE]
  sigma <- x$sigma[filas]
  indicadores <- x$X[, c("zona_urbana", "zona_aldea"), drop = FALSE]
  stopifnot(all(indicadores %in% c(0, 1)), all(rowSums(indicadores) <= 1))
  zona <- ifelse(indicadores[, 1] == 1, "Urbana",
                 ifelse(indicadores[, 2] == 1, "Aldea", "Rural"))
  focal <- zona == zona_focal & (is.na(region_focal) || x$codigo == region_focal)
  por_zona <- array(0, c(n_replicas, length(zonas), length(cobertura)),
                    dimnames = list(NULL, zonas, names(cobertura)))

  set.seed(semilla)
  for (inicio in seq.int(1L, length(x$V), by = 512L)) {
    j <- inicio:min(inicio + 511L, length(x$V))
    # 1. Predictor por réplica y segmento. El índice de comuna reemplaza Z %*% gamma.
    eta <- beta %*% t(x$X[j, , drop = FALSE]) + gamma[, x$grupo_comuna[j], drop = FALSE]
    # 2. Densidad simulada y limitada por el techo del modelo.
    z <- matrix(rnorm(length(eta)), nrow = n_replicas)
    densidad <- pmin(exp(eta + z * sigma), x$delta_max)
    u <- matrix(runif(length(eta)), nrow = n_replicas)
    coberturas <- list(base = rep(1, length(j)), uniforme = rep(cobertura[["uniforme"]], length(j)),
      focalizada = ifelse(focal[j], cobertura[["focalizada"]], 1))
    # La base usa V_tot. Los segmentos que mantienen exposición reutilizan su conteo base.
    lambda_base <- pmin(pmax(sweep(densidad, 2L, x$V[j], "*"), 1), 1e7)
    y_base <- matrix(qpois(as.vector(u), as.vector(lambda_base)), nrow = n_replicas)
    for (e in names(coberturas)) {
      # 3. Exposición del escenario.
      viviendas <- if (e == "base") x$V[j] else round(x$V[j] / coberturas[[e]])
      cambia <- viviendas != x$V[j]
      y <- y_base
      if (any(cambia)) {
        # 4. Intensidad con los límites de Stan y conteo Poisson con el mismo U.
        lambda <- pmin(pmax(sweep(densidad[, cambia, drop = FALSE], 2L,
                                  viviendas[cambia], "*"), 1), 1e7)
        y[, cambia] <- qpois(as.vector(u[, cambia, drop = FALSE]), as.vector(lambda))
      }
      # 5. Sumar dentro de cada réplica en la zona a la que pertenece el segmento.
      for (g in zonas) {
        por_zona[, g, e] <- por_zona[, g, e] + rowSums(y[, zona[j] == g, drop = FALSE])
      }
    }
  }
  totales <- vapply(names(cobertura), function(e) rowSums(por_zona[, , e]), numeric(n_replicas))
  stopifnot(all(is.finite(totales)), all(totales >= 0))
  list(totales = totales, por_zona = por_zona, segmentos = length(x$V),
       segmentos_zona = table(factor(zona, levels = zonas)))
}

#' Convertir un resumen predictivo en una fila
#'
#' @param x Vector con una suma o diferencia para cada réplica.
#' @return Tabla de una fila con media y percentiles 2,5, 50 y 97,5.
tabla_resumen <- function(x) as.data.frame(t(resumir_replicas(x)))

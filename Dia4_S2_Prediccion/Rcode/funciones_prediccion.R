# Funciones compartidas de D4S2. Cargar con source().

#' Construir llaves de segmento
#'
#' @param d Tabla con las columnas edid y fid.
#' @return Vector de llaves en el orden de las filas.
key <- function(d) paste(d$edid, d$fid, sep = "\r")

#' Alinear una tabla de coberturas
#'
#' @param tabla Tabla con edid, fid y cobertura, una fila válida por segmento.
#' @param destino Tabla de segmentos en el orden que usa la predicción.
#' @return Vector numérico de coberturas en el orden del destino.
alinear_cobertura <- function(tabla, destino) {
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

#' Resumir totales por réplica de Magallanes
#'
#' @param x Vector con un conteo total por réplica.
#' @return Vector nombrado con media, mediana, percentil 2,5 y percentil 97,5.
#' @examples
#' summarize_total(c(150, 160, 175))
summarize_total <- function(x) {
  c(media = mean(x), mediana = median(x),
    inferior = unname(quantile(x, 0.025)), superior = unname(quantile(x, 0.975)))
}

#' Formatear las etiquetas de un eje
#'
#' @param x Valores numéricos que se mostrarán en el eje.
#' @return Etiquetas con punto de miles y coma decimal.
spanish_labels <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)

#' Leer los totales originales de Stan para una región
#'
#' @param codigo Código regional de dos dígitos.
#' @return Lista con las réplicas ordenadas por cadena e iteración y una fila
#'   de metadatos del destino, incluidos los segmentos con exposición cero.
leer_replica_regional <- function(codigo) {
  x <- leer_region(codigo)
  stopifnot(!anyNA(x$cadena), !anyNA(x$iteracion),
    !anyDuplicated(data.frame(cadena = x$cadena, iteracion = x$iteracion)),
    all(is.finite(x$total_stan)), all(x$total_stan >= 0),
    all(x$total_stan == floor(x$total_stan)))
  list(replicas = x$total_stan[order(x$cadena, x$iteracion)],
    destino = data.frame(codigo = codigo, region = x$region,
      segmentos = length(x$V), exposicion_cero = sum(x$V == 0),
      observado = x$observado))
}

#' Resumir la distribución de un total
#'
#' @param x Vector con un conteo total por réplica.
#' @return Vector nombrado con media, percentil 2,5, mediana y percentil 97,5.
#' @examples
#' resumir_total(c(150, 160, 175))
resumir_total <- function(x) {
  q <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  c(media = mean(x), inferior = q[1], mediana = q[2], superior = q[3])
}

#' Calcular los percentiles de un total por zona
#'
#' @param x Vector con un conteo total por réplica.
#' @return Vector con los percentiles 2,5, 50 y 97,5.
#' @examples
#' cuantiles_zona(c(150, 160, 175))
cuantiles_zona <- function(x) unname(quantile(x, c(0.025, 0.5, 0.975)))

#' Regenerar las predicciones por zona de una región
#'
#' @param codigo Código regional de dos dígitos.
#' @param pos Posiciones de las réplicas en el emparejamiento nacional.
#' @param referencia Tabla local con segmentos, censados y medias originales por región y zona.
#' @param ajuste Agregado nacional con las posiciones y los totales de Stan.
#' @param categorias Nombres de zona usados como columnas de salida.
#' @return Lista con la tabla regional por zona y una matriz de réplicas por zona.
#' @details Fija la semilla por región para repetir los conteos nuevos.
simular_zonas <- function(codigo, pos, referencia, ajuste, categorias) {
  x <- leer_region(codigo)
  ref <- referencia[referencia$codigo == codigo, , drop = FALSE]
  stopifnot(nrow(ref) == length(categorias), !anyDuplicated(ref$zona),
    setequal(ref$zona, categorias), sum(ref$censado) == x$observado,
    identical(as.numeric(x$total_stan[pos]), as.numeric(ajuste$replicas_regionales[, codigo])))
  indicadores <- x$X[, c("zona_urbana", "zona_aldea"), drop = FALSE]
  stopifnot(all(indicadores %in% c(0, 1)), all(rowSums(indicadores) <= 1))
  zona <- ifelse(indicadores[, 1] == 1, "Urbana",
                 ifelse(indicadores[, 2] == 1, "Aldea", "Rural"))
  beta <- x$beta[pos, , drop = FALSE]
  gamma <- x$gamma[pos, , drop = FALSE]
  sigma <- x$sigma[pos]
  m <- length(pos)
  totales <- matrix(0, m, length(categorias), dimnames = list(NULL, categorias))
  set.seed(2700L + as.integer(codigo))
  # Procesar segmentos por bloques para limitar memoria sin partir las réplicas.
  for (inicio in seq.int(1L, nrow(x$X), by = 512L)) {
    j <- inicio:min(inicio + 511L, nrow(x$X))
    eta <- beta %*% t(x$X[j, , drop = FALSE]) + gamma[, x$grupo_comuna[j], drop = FALSE]
    densidad <- pmin(exp(eta + matrix(rnorm(length(eta)), m) * sigma), x$delta_max)
    # Aplicar el techo de densidad y los límites de intensidad del programa Stan.
    lambda <- pmin(pmax(sweep(densidad, 2, x$V[j], "*"), 1), 1e7)
    y <- matrix(rpois(length(lambda), lambda), m)
    for (z in categorias) {
      totales[, z] <- totales[, z] + rowSums(y[, zona[j] == z, drop = FALSE])
    }
  }
  # Sumar las medias originales. Los cuantiles provienen de las nuevas réplicas.
  filas <- lapply(categorias, function(z) {
    i <- zona == z
    q <- cuantiles_zona(totales[, z])
    fila <- ref[ref$zona == z, , drop = FALSE]
    stopifnot(sum(i) == fila$segmentos)
    data.frame(codigo = codigo, zona = z, segmentos = fila$segmentos,
      censado = fila$censado, media_original = fila$media_original,
      inferior = q[1], mediana = q[2], superior = q[3])
  })
  list(tabla = do.call(rbind, filas), replicas = totales)
}

#' Cargar el paquete docente de una región
#'
#' @param codigo Código regional de dos dígitos entre "01" y "16".
#' @return Lista con el diseño, la exposición, las cadenas y los totales originales.
#' @details Busca el archivo en data/regiones desde la carpeta de la práctica.
leer_region <- function(codigo) {
  stopifnot(length(codigo) == 1, codigo %in% sprintf("%02d", 1:16))
  x <- readRDS(file.path("data/regiones", paste0("region_", codigo, ".rds")))
  stopifnot(identical(x$codigo, codigo), nrow(x$X) == length(x$V),
    length(x$grupo_comuna) == length(x$V),
    all(x$grupo_comuna %in% seq_len(ncol(x$gamma))),
    ncol(x$X) == ncol(x$beta), nrow(x$beta) == length(x$sigma),
    nrow(x$gamma) == length(x$sigma), length(x$total_stan) == length(x$sigma),
    all(x$V >= 0), all(x$sigma > 0))
  x
}

#' Simular conteos para un bloque de segmentos
#'
#' @param viviendas Vector de exposiciones del bloque.
#' @param densidad Matriz de densidades simuladas.
#' @param u Matriz de uniformes compartidos entre escenarios.
#' @param m Número de réplicas.
#' @return Matriz de conteos Poisson con una columna por segmento.
simular_conteos_region <- function(viviendas, densidad, u, m) {
  lambda <- pmin(pmax(sweep(densidad, 2, viviendas, "*"), 1), 1e7)
  matrix(qpois(u, lambda), m)
}

#' Predecir una región con dos escenarios de cobertura
#'
#' @param x Paquete validado por leer_region().
#' @param n_replicas Número de réplicas guardadas que se usarán.
#' @param cobertura Cobertura supuesta en (0, 1]. Use 1 para conservar la exposición base.
#' @param semilla Semilla para repetir la simulación.
#' @return Matriz con totales por réplica para la base y el escenario. El atributo por_zona contiene las réplicas desglosadas.
#' @details Fija la semilla de R. Reutiliza las densidades y los uniformes para comparar ambas exposiciones.
predecir_region <- function(x, n_replicas = 500L, cobertura = 0.95, semilla = 2625L) {
  stopifnot(length(n_replicas) == 1, is.finite(n_replicas),
    n_replicas == floor(n_replicas), n_replicas >= 2, n_replicas <= length(x$sigma),
    length(cobertura) == 1, is.finite(cobertura), cobertura > 0, cobertura <= 1)
  # Seleccionar filas espaciadas de las cadenas guardadas, manteniendo juntos beta, gamma y sigma.
  filas <- unique(round(seq(1, length(x$sigma), length.out = n_replicas)))
  m <- length(filas)
  beta <- x$beta[filas, , drop = FALSE]
  gamma <- x$gamma[filas, , drop = FALSE]
  sigma <- x$sigma[filas]
  zonas <- c("Rural", "Aldea", "Urbana")
  escenario <- c("Base en R", "Escenario en R")
  indicadores <- x$X[, c("zona_urbana", "zona_aldea"), drop = FALSE]
  stopifnot(all(indicadores %in% c(0, 1)), all(rowSums(indicadores) <= 1))
  zona <- ifelse(indicadores[, 1] == 1, "Urbana",
                 ifelse(indicadores[, 2] == 1, "Aldea", "Rural"))
  por_zona <- array(0, c(m, length(zonas), length(escenario)),
                    dimnames = list(NULL, zonas, escenario))
  set.seed(semilla)
  # Cada fila usa el mismo conjunto de parámetros para todos los segmentos.
  for (inicio in seq.int(1L, length(x$V), by = 512L)) {
    j <- inicio:min(inicio + 511L, length(x$V))
    eta <- beta %*% t(x$X[j, , drop = FALSE]) +
      gamma[, x$grupo_comuna[j], drop = FALSE]
    densidad <- pmin(exp(eta + matrix(rnorm(length(eta)), m) * sigma), x$delta_max)
    # Usar los mismos uniformes para que la comparación aísle el cambio de exposición.
    u <- matrix(runif(length(eta)), m)
    # Escalar las viviendas por la cobertura supuesta y volver a simular los conteos.
    base <- simular_conteos_region(x$V[j], densidad, u, m)
    alternativa <- simular_conteos_region(round(x$V[j] / cobertura), densidad, u, m)
    for (z in zonas) {
      i <- zona[j] == z
      por_zona[, z, 1] <- por_zona[, z, 1] + rowSums(base[, i, drop = FALSE])
      por_zona[, z, 2] <- por_zona[, z, 2] + rowSums(alternativa[, i, drop = FALSE])
    }
  }
  totales <- vapply(escenario, function(e) rowSums(por_zona[, , e]), numeric(m))
  stopifnot(all(is.finite(totales)), all(totales[, 2] >= totales[, 1]))
  attr(totales, "por_zona") <- por_zona
  totales
}

## ---- prediccion-funcion ----
#' Predecir los segmentos de Magallanes
#'
#' @param cobertura Vector de coberturas en (0, 1], una por segmento.
#' @param destino Tabla ordenada de segmentos con la exposición V_tot.
#' @param density Matriz de densidades simuladas, con réplicas en filas.
#' @param uniforms Matriz de uniformes con las mismas dimensiones que density.
#' @return Lista con conteos, intensidades esperadas, viviendas y coberturas.
#' @examples
#' predecir(1, data.frame(V_tot = 50), matrix(3), matrix(0.5))
predecir <- function(cobertura, destino, density, uniforms) {
  stopifnot(length(cobertura) == nrow(destino), all(is.finite(cobertura)),
            all(cobertura > 0 & cobertura <= 1))
  viviendas <- round(destino$V_tot / cobertura)
  lambda <- pmin(pmax(sweep(density, 2, viviendas, "*"), 1), 1e7)
  # qpois transforma los mismos uniformes en conteos Poisson por escenario.
  y <- matrix(qpois(as.vector(uniforms), as.vector(lambda)), nrow(lambda))
  stopifnot(all(is.finite(y)), all(y >= 0))
  list(y = y, lambda = lambda, viviendas = viviendas, cobertura = cobertura)
}

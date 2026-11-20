# D4S2 | Agregación nacional de las predicciones originales de Stan
# Ejecución: desde Dia4_S2_Prediccion, ejecutar este script primero.
# Insumos: data/regiones/region_XX.rds, con las réplicas originales de Stan.
# Salidas: prediccion_16_regiones.csv, prediccion_chile.csv y
# prediccion_chile.rds, con totales por réplica e índices de emparejamiento.
# El destino del modelo es un subconjunto censal y las cifras son preliminares.

# Ejecutar desde Dia4_S2_Prediccion. Usa las predicciones regionales originales.
source("Rcode/funciones_prediccion.R",
       local = TRUE,
       encoding = "UTF-8")

codigos <- sprintf("%02d", 1:16)

## ---- nacional-leer ----
# Leer los totales de Stan y recuperar su orden de cadena e iteración.
insumos <- setNames(lapply(codigos, leer_replica_regional), codigos)
replicas_regionales <- lapply(insumos, function(x) x$replicas)

## ---- nacional-agregar ----
# Los modelos se ajustaron por región. Se combina una réplica aleatoria de cada una.
# La permutación conserva todas las réplicas marginales, sin alinear cadenas por número.
n_replicas <- unique(lengths(replicas_regionales))
stopifnot(length(n_replicas) == 1, n_replicas > 1L)
semilla <- 2614L
set.seed(semilla)
# Cada columna permuta las réplicas de una región de forma independiente.
indice <- vapply(codigos, function(codigo)
  sample.int(n_replicas), integer(n_replicas))
conjunto <- vapply(seq_along(codigos), function(j) {
  replicas_regionales[[j]][indice[, j]]
}, numeric(n_replicas))
colnames(conjunto) <- codigos
# Paso 07 del proyecto: recuento fijo de personas censadas sin segmento.
# Fuente del recuento: src/07_estimacion_edificaciones_geoespaciales.R del proyecto CHL2026.
sin_segmento <- 254224L
total_segmentos <- rowSums(conjunto)
total_chile <- total_segmentos + sin_segmento


## ---- nacional-resumir ----
# Calcular percentiles después de sumar dentro de cada réplica conjunta.
por_region <- as.data.frame(t(vapply(
  replicas_regionales, resumir_total, numeric(4)
)))
stopifnot(abs(mean(total_segmentos) - sum(por_region$media)) < 1e-7)
destinos <- do.call(rbind, lapply(insumos, function(x)
  x$destino))
por_region <- cbind(destinos, por_region)
# Mantener por separado el universo con segmento y el agregado con el recuento fijo.
pais_segmentos <- cbind(
  data.frame(
    universo = "Segmentos destino de las 16 regiones",
    segmentos = sum(por_region$segmentos),
    exposicion_cero = sum(por_region$exposicion_cero),
    observado = sum(por_region$observado)
  ),
  as.data.frame(t(resumir_total(total_segmentos)))
)

pais <- cbind(
  data.frame(
    universo = "Segmentos destino y censo sin segmento",
    segmentos = pais_segmentos$segmentos,
    exposicion_cero = pais_segmentos$exposicion_cero,
    personas_sin_segmento = sin_segmento,
    observado = pais_segmentos$observado + sin_segmento
  ),
  as.data.frame(t(resumir_total(total_chile)))
)
print(pais)

dir.create("resultados", showWarnings = FALSE)
write.csv(
  por_region,
  "resultados/prediccion_16_regiones.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  pais,
  "resultados/prediccion_chile.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
saveRDS(
  list(
    pais = pais,
    pais_segmentos = pais_segmentos,
    regiones = por_region,
    total_chile = total_chile,
    total_segmentos = total_segmentos,
    replicas_regionales = conjunto,
    indice_regional = indice,
    semilla = semilla
  ),
  "resultados/prediccion_chile.rds"
)

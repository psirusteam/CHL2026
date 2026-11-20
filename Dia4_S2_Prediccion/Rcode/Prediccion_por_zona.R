# D4S2 | Predicciones para las zonas rural, aldea y urbana
# Ejecución: desde Dia4_S2_Prediccion, después de Prediccion_nacional.R.
# Insumos: prediccion_chile.rds, data/regiones y data/referencia_region_zona.rds.
# Salidas: prediccion_por_zona.csv (país), prediccion_region_zona.csv
# (región y zona) y prediccion_por_zona.rds (también contiene réplicas).
# Las medias de Stan son originales. Las medianas e intervalos se obtienen de
# conteos nuevos en R. El censo sin segmento se agrega solo al total nacional.

# Ejecutar después de Prediccion_nacional.R, desde Dia4_S2_Prediccion.
source("Rcode/funciones_prediccion.R", local = TRUE, encoding = "UTF-8")

ajuste <- readRDS("resultados/prediccion_chile.rds")
categorias <- c("Rural", "Aldea", "Urbana")
referencia <- readRDS("data/referencia_region_zona.rds")
codigos <- names(ajuste$replicas_regionales[1, ])
resultados <- lapply(seq_along(codigos), function(i) {
  cat("Predicción por zona, región", codigos[i], "\n")
  flush.console()
  simular_zonas(codigos[i], ajuste$indice_regional[, i], referencia, ajuste, categorias)
})
regional <- do.call(rbind, lapply(resultados, function(x) x$tabla))
regional$region <- ajuste$regiones$region[match(regional$codigo, ajuste$regiones$codigo)]
# Cada fila combina las mismas posiciones emparejadas en el agregado nacional.
sumas <- Reduce("+", lapply(resultados, function(x) x$replicas))
nacional <- do.call(rbind, lapply(categorias, function(z) {
  s <- regional[regional$zona == z, ]
  q <- cuantiles_zona(sumas[, z])
  data.frame(zona = z, segmentos = sum(s$segmentos), censado = sum(s$censado),
    media_original = sum(s$media_original), inferior = q[1], mediana = q[2], superior = q[3])
}))
stopifnot(sum(nacional$segmentos) == ajuste$pais_segmentos$segmentos,
  sum(nacional$censado) == ajuste$pais_segmentos$observado,
  isTRUE(all.equal(sum(nacional$media_original), ajuste$pais_segmentos$media)))
print(nacional)
write.csv(nacional, "resultados/prediccion_por_zona.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(regional, "resultados/prediccion_region_zona.csv", row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(nacional = nacional, regiones = regional, replicas = sumas),
        "resultados/prediccion_por_zona.rds")

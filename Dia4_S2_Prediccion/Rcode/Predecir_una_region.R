# D4S2 | Práctica reproducible de predicción regional con R base
# Ejecución: abrir Practica_D4S2.Rproj o fijar la carpeta de la práctica.
# Elegir la región con Sys.setenv(CHL2026_REGION = "05"), por ejemplo, y ejecutar
# source("Rcode/Predecir_una_region.R", encoding = "UTF-8").
# Insumos: data/regiones/region_XX.rds, con covariables, exposición y cadenas.
# Salida: tablas en la consola del total regional y de cada zona, con media,
# mediana e intervalo predictivo del 95 %. No requiere red ni compilar Stan.
# La cobertura es un supuesto de la práctica. Los resultados son preliminares.

# Abrir el proyecto o fijar el directorio de trabajo en esta carpeta.
# Cambiar el código entre "01" y "16". Se usa R base y las cadenas ya ajustadas.
source("Rcode/funciones_prediccion.R", local = TRUE, encoding = "UTF-8")

codigo_region <- Sys.getenv("CHL2026_REGION", unset = "12")
n_replicas <- 500L
cobertura <- 0.95
semilla <- 2625L

# Recibir un total por réplica y resumir su distribución predictiva.
# Cargar una región y verificar correspondencias entre segmentos y coeficientes.
region <- leer_region(codigo_region)
totales_r <- predecir_region(region, n_replicas, cobertura, semilla)
# Comparar la referencia original de Stan con las dos simulaciones nuevas en R.
tabla <- rbind("Base original de Stan" = resumir_total(region$total_stan),
              t(apply(totales_r, 2, resumir_total)))
cat(region$region, "-", length(region$V), "segmentos del destino\n")
cat("Personas censadas en el mismo destino:", region$observado, "\n")
print(round(tabla))
cat("\nResultados simulados en R por zona:\n")
por_zona <- attr(totales_r, "por_zona")
for (escenario in dimnames(por_zona)[[3]]) {
  cat(escenario, "\n")
  print(round(t(apply(por_zona[, , escenario], 2, resumir_total))))
}
# La base en R usa una simulación nueva. Comparar resúmenes, no valores fila a fila.
# Con las alertas de D4S1, los resultados son preliminares.

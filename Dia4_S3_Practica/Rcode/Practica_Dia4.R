# Práctica D4S3. Abrir el proyecto o fijar esta carpeta en RStudio.
# Ejecutar los bloques en orden. Se usan cadenas ya estimadas y R base.
# La ejecución completa escribe tablas y réplicas en resultados/.

## ---- ejercicio-configuracion ----
source("Rcode/funciones_practica.R", encoding = "UTF-8")
codigos <- sprintf("%02d", 1:16)
n_replicas <- 1000L
cobertura <- c(base = 1, uniforme = 0.95, focalizada = 0.90)
zona_focal <- "Rural"           # También puede ser "Aldea" o "Urbana".
region_focal <- NA_character_   # Todas las regiones. Usar "12" solo para Magallanes.
semilla <- 2631L
semilla_pais <- 2614L

## ---- ejercicio-manual ----
# Una réplica de un segmento. Mantener z y u al cambiar la cobertura.
eta_manual <- log(3)
sigma_manual <- 0.2
z_manual <- 0.5
u_manual <- 0.5
densidad_manual <- min(exp(eta_manual + sigma_manual * z_manual), 4.35)
c_manual <- c(1, 0.90)
v_manual <- round(90 / c_manual)
lambda_manual <- pmin(pmax(densidad_manual * v_manual, 1), 1e7)
print(data.frame(cobertura = c_manual, viviendas = v_manual,
  densidad = densidad_manual, intensidad = lambda_manual,
  personas = qpois(u_manual, lambda_manual)))

# Cinco réplicas de dos segmentos. Ordenar las sumas y calcular percentiles.
y_manual <- cbind(c(63, 68, 72, 75, 82), c(42, 44, 47, 46, 47))
t_manual <- sort(rowSums(y_manual))
p_manual <- c(0.025, 0.50, 0.975)
h_manual <- 1 + (length(t_manual) - 1) * p_manual
j_manual <- floor(h_manual)
a_manual <- h_manual - j_manual
q_manual <- (1 - a_manual) * t_manual[j_manual] + a_manual * t_manual[j_manual + 1L]
print(data.frame(percentil = p_manual, posicion = h_manual, cuantil = q_manual))
stopifnot(isTRUE(all.equal(q_manual, unname(quantile(t_manual, p_manual, type = 7)))))

## ---- ejercicio-referencia ----
insumos <- setNames(lapply(codigos, leer_insumo_regional), codigos)
n_stan <- unique(vapply(insumos, function(x) length(x$total_stan), integer(1)))
stopifnot(length(n_stan) == 1L)
indices_stan <- emparejar_regiones(n_stan, codigos, semilla_pais)
totales_stan <- vapply(codigos, function(codigo) {
  x <- insumos[[codigo]]
  x$total_stan[order(x$cadena, x$iteracion)][indices_stan[, codigo]]
}, numeric(n_stan))
censado_destino <- sum(vapply(insumos, function(x) as.numeric(x$observado), numeric(1)))
# Paso 07 del proyecto CHL2026. Este conteo fijo se suma solo a Chile.
personas_sin_segmento <- 254224L
referencia <- cbind(data.frame(censado_destino, personas_sin_segmento),
  as.data.frame(t(resumir_replicas(rowSums(totales_stan) + personas_sin_segmento))))
print(referencia)

## ---- ejercicio-region ----
codigo_practica <- "12"
ejemplo <- simular_coberturas(insumos[[codigo_practica]], n_replicas, cobertura,
  zona_focal, region_focal, semilla + as.integer(codigo_practica))
print(round(t(apply(ejemplo$totales, 2, resumir_replicas))))

## ---- ejercicio-nacional ----
# El mismo emparejamiento regional se usa para los tres escenarios.
idx <- emparejar_regiones(n_replicas, codigos, semilla_pais)
resultados <- lapply(codigos, function(codigo) {
  message("Procesando región ", codigo)
  if (codigo == codigo_practica) return(ejemplo)
  simular_coberturas(insumos[[codigo]], n_replicas, cobertura, zona_focal,
                    region_focal, semilla + as.integer(codigo))
})
names(resultados) <- codigos
por_zona_chile <- Reduce("+", lapply(codigos, function(codigo) {
  resultados[[codigo]]$por_zona[idx[, codigo], , , drop = FALSE]
}))
sin_agregar <- apply(por_zona_chile, c(1L, 3L), sum)
total_chile <- sin_agregar + personas_sin_segmento
print(round(t(apply(total_chile, 2L, resumir_replicas))))

## ---- ejercicio-guardar ----
escenarios <- names(cobertura)
zonas_destino <- dimnames(por_zona_chile)[[2L]]
regiones <- do.call(rbind, lapply(codigos, function(codigo) {
  x <- insumos[[codigo]]
  do.call(rbind, lapply(escenarios, function(e) {
    cbind(data.frame(codigo = codigo, region = x$region, escenario = e,
      segmentos = length(x$V), censado = x$observado),
      tabla_resumen(resultados[[codigo]]$totales[, e]))
  }))
}))
zonas <- do.call(rbind, lapply(escenarios, function(e) {
  do.call(rbind, lapply(zonas_destino, function(zona) {
    cbind(data.frame(zona, escenario = e), tabla_resumen(por_zona_chile[, zona, e]))
  }))
}))
pais <- do.call(rbind, lapply(escenarios, function(e) {
  cbind(data.frame(escenario = e, censado_destino, personas_sin_segmento),
        tabla_resumen(total_chile[, e]))
}))
# Restar las réplicas emparejadas y calcular los percentiles de las diferencias.
contraste <- do.call(rbind, lapply(escenarios[-1L], function(e) {
  cbind(data.frame(escenario = e), tabla_resumen(total_chile[, e] - total_chile[, "base"]))
}))
diag <- read.csv("data/diagnosticos_regionales.csv", fileEncoding = "UTF-8",
                 colClasses = c(codigo = "character"))
stopifnot(setequal(diag$codigo, codigos), !anyDuplicated(diag$codigo))
alertas <- sum(diag$parametros_con_alerta > 0)
dir.create("resultados", showWarnings = FALSE)
write.csv(regiones, "resultados/escenarios_regionales.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(zonas, "resultados/escenarios_por_zona.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pais, "resultados/escenarios_chile.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(contraste, "resultados/diferencias_chile.csv", row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(regiones = regiones, zonas = zonas, pais = pais, contraste = contraste,
    referencia = referencia, alertas = alertas, cobertura = cobertura, n_replicas = n_replicas,
    zona_focal = zona_focal, region_focal = region_focal, semilla = semilla, semilla_pais = semilla_pais,
    codigo_practica = codigo_practica, replicaciones_region = resultados,
    replicas_pais = total_chile, replicas_zona = por_zona_chile, indices_region = idx),
  "resultados/resumen_practica.rds")
capture.output(sessionInfo(), file = "resultados/sessionInfo_practica.txt")

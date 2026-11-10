# Modelo final de Chile: referencia común del taller

Los paquetes `region_11.rds` y `region_12.rds` exportan los ajustes implementados de Aysén y Magallanes desde `output/04_seleccion_covariables_generalizado/cut_zona_fija/ajuste`. No se ejecutó otro muestreo para producirlos. El Stan tiene firma MD5 `2babb95f5c0bfb2565ebe63fa9a60435`.

| Región | Aptos | Destino | Columnas X | Comunas | Covariables |
|---|---:|---:|---:|---:|---:|
| Aysén | 1455 | 2085 | 43 | 10 | 40 |
| Magallanes | 2288 | 3044 | 43 | 10 | 40 |

Se preservan todas las iteraciones posteriores en su orden de cadena (500 por cadena, cuatro cadenas), parámetros, controles NUTS, cantidades generadas originales, matrices, llaves, exposición y configuración. Las covariables `_std` se toman de X_tot sin reestandarizar. El orden de zona es intercepto, urbana, aldea y las 40 covariables de `covariables_modelo.csv`. Los niveles Z conservan sus nombres `cut_NNNNN` y el diccionario numérico de comunas por separado.

## Cómo continuar

1. Ejecutar `Rcode/Modelo_Chile.R` desde D3S3 para Aysén y `Rcode/Practica_Dia3.R` desde D3S4 para Magallanes. Por defecto se usa el ajuste original. `reestimar = TRUE` produce una corrida opcional separada, que no sustituye la referencia.
2. Ejecutar `Rcode/Diagnostico.R` desde D4S1. `Preparar_cadenas_docente.R` vuelve a importar el paquete original, sin muestrear.
3. Ejecutar `Rcode/Prediccion.R` desde D4S2 y `Rcode/Practica_Dia4.R` desde D4S3. Ambos cotejan la procedencia, columnas, grupos y llaves con D4S1.
4. Regenerar las presentaciones con Quarto desde sus carpetas.

Los escenarios adicionales usan `round(V_tot_final / c_adicional)`. Con c=1 se conserva la exposición del modelo final. Las reglas de techo de densidad y piso/tope de intensidad son exactamente las implementadas, incluidos los casos con exposición cero. Replicar esas reglas no valida sus supuestos.

## Estado del ajuste y lectura de resultados

Los ajustes originales tienen R-hat máximos próximos a 3,36 (Aysén) y 3,30 (Magallanes), ESS bulk mínimos próximos a 4,5 y rango X de 42 para 43 columnas. Las alertas se conservan. Las tablas y escenarios son ejercicios explícitos de cálculo pendientes de validación; no se deben comunicar como estimaciones poblacionales validadas. La nueva actividad de colinealidad y los ejercicios de diagnóstico discuten qué revisar antes de una nueva especificación. No se eliminaron variables ni se cambiaron las previas para ocultar este comportamiento.

## Regeneración desde el proyecto

Desde `Taller/Revisiones/v02_2026-10-21`, ejecutar `Rscript Modelo_final/Preparar_insumos.R`. Requiere acceso al proyecto y memoria suficiente para leer un stanfit regional. El taller distribuido usa los paquetes compactos y no necesita cargar los stanfit de varios GB. Los resultados originales del proyecto no se modifican.

`Rscript Modelo_final/Validar_insumos.R` verifica equivalencia de matrices, exposición, llaves, parámetros y totales originales, y coherencia entre las sesiones. La trazabilidad conserva las firmas de los insumos fuente en `configuration$input_hashes` y las configuraciones originales dentro de cada paquete.

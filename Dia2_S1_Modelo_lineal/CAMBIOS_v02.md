# D2S1 · Introducción a los modelos lineales con métodos bayesianos

Versión de trabajo v02 · Revisión del 23 de octubre de 2026.

## Programa y recorrido

Martes 24 de noviembre de 2026, 09:30–11:00. Docente: Andrés Gutiérrez. Se mantienen regresión de log(personas por vivienda), covariables, comparación con MCO, ajuste Stan, sensibilidad, R² bayesiano, chequeo predictivo y predicción. Se incorpora la conexión con D1S3 y D1S4 y un recorrido orientativo de 90 minutos; D2S2 retoma la práctica.

## Cambios conceptuales

- Intercepto interpretado en la zona rural y en las medias de las covariables continuas, sin equipararlo con el promedio nacional de la respuesta.
- Notación p coherente como número de columnas de diseño, incluido el intercepto.
- Efectos multiplicativos presentados sobre medianas condicionales; distinción entre una desviación estándar continua y un cambio de categoría.
- Estandarización: invariancia de MCO bajo cambios invertibles, pero necesidad de transformar también las previas para preservar un modelo bayesiano equivalente.
- Elegibilidad: se conserva el filtro operativo apto_modelo. Su selección no convierte automáticamente la verosimilitud Normal en truncada, y la predicción fuera del subconjunto requiere evaluar transportabilidad.
- Regularización: sigma fija y penalización también del intercepto en este ejemplo. El encogimiento actúa en direcciones del diseño; no necesariamente reduce cada coeficiente individual.
- Separación entre condicional de beta con sigma fijada en su estimador y posterior marginal de Stan, que integra sigma. La cercanía de estimaciones no prueba por sí sola la corrección del muestreo.
- Se retiran afirmaciones automáticas de convergencia, adecuación de previas y atribución causal de problemas del modelo de Chile a colinealidad.
- R² y chequeos predictivos presentados como evaluaciones parciales del ajuste, sin asegurar que una densidad global favorable valida los subgrupos.
- Aproximación delta del ruido Poisson con sus condiciones; la comparación descriptiva no constituye una descomposición formal de la varianza.
- Predicción de densidad diferenciada de población corregida por cobertura; se mantiene la escala de covariables del ajuste.

## Cambios de código

- Comprobación de densidades positivas y finitas, diseño finito, n > p y rango completo.
- Asignación de nombres de coeficientes por sus identificadores Stan.
- Media posterior del efecto porcentual calculada draw a draw, en lugar de transformar la media de beta.
- Condicional evaluada en sigma fija identificada como condicional_sigma_fija.
- Semillas distintas y calentamiento explícito para el ajuste base y el de sensibilidad.
- ESS_bulk, ESS_tail, error Monte Carlo, divergencias y profundidad del árbol para ambos ajustes; el control de profundidad usa el límite predeterminado de 10 de estas ejecuciones.
- Registro de sessionInfo() y verificación de la matriz de diseño utilizada al predecir.

## Trazabilidad y conservación

La elegibilidad se contrastó con funciones/fun_modelo_poblacion.R; la estandarización con src/02_estandarizacion_covariables.R y funciones/fun_estandarizacion.R. Se conservan los datos originales y una sola versión de trabajo de la sesión. El material inicial permanece intacto en Taller/Dia2_S1_Modelo_lineal.

Los diagnósticos siguen las referencias oficiales de Stan y posterior consultadas durante la revisión de D1S3. Los detalles de ejecución y compilación quedan en resultados_log.txt y render_log.txt. El resultado final se añade al terminar.

### Resultado de las verificaciones

- Ambos ajustes se ejecutaron; se revisaron 14 parámetros de cada uno. R-hat máximo 1,00354, ESS_bulk mínimo superior a 1700 y ESS_tail mínimo superior a 2200. No se registraron divergencias ni alcances de profundidad máxima.
- Se verificó que el R² posterior estuviera entre cero y uno y que el código final pudiera analizarse sin errores de sintaxis.
- La curva exploratoria LOESS produjo avisos por vecindades de ancho cero y valores repetidos. Se sustituyó por una tendencia lineal descriptiva y se ejecutó nuevamente ese bloque; el registro específico queda en verificacion_final_log.txt. Los ajustes bayesianos no cambiaron.
- El intercepto no se reporta como efecto porcentual: su celda queda sin valor en esa columna. El código y la tabla generada se actualizaron de forma coherente.
- Se separaron código y resultados de la comparación para reducir la carga por lámina.
- La ejecución registró avisos de configuración regional C.UTF-8 que no impidieron producir los resultados. La revisión visual completa de las diapositivas queda pendiente.

- Compilación final: HTML generado sin avisos de estructura de Quarto. Enlaces del README verificados.


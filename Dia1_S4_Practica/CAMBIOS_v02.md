# D1S4 · Práctica computacional en RStan

Versión de trabajo v02 · Revisión del 23 de octubre de 2026.

## Correspondencia con el programa

Lunes 23 de noviembre, 16:30–18:00. Docente: Stalyn Guerrero. La práctica replica los fundamentos uniparamétricos y multiparamétricos con datos de Magallanes y compara con Aysén, según el recorrido de D1S2 y D1S3. Identidad institucional y logos alineados con las sesiones revisadas.

## Organización didáctica

- Ejercicios 1–3 como recorrido principal: Bernoulli–Beta, Poisson–Gamma y Normal por zona.
- Ejercicios 4 y 5 conservados como extensiones: Multinomial–Dirichlet y desafío Metropolis.
- Orientación temporal para los 90 minutos y producto concreto por pareja: estimaciones e intervalos, diagnóstico, chequeo predictivo y conclusión.
- Referencias a las sesiones previas corregidas: los modelos multiparamétricos corresponden también a la sesión de la tarde.

## Correcciones estadísticas y de código

1. La fila de Stan del ejercicio 1 ahora calcula la media posterior; antes usaba el cuantil 0,5 bajo una columna de medias. Se añadió el ajuste Stan de sensibilidad Beta(1,99), además de su solución analítica.
2. Se distingue la proporción descriptiva de un conjunto observado de la probabilidad inferida bajo Bernoulli/Multinomial para segmentos comparables.
3. Se evita equiparar la presencia de viviendas colectivas con el indicador operativo excluido_colectivo, y se retira la afirmación de que este ejercicio determina cómo predecir población colectiva.
4. La independencia posterior entre regiones requiere verosimilitudes y previas independientes, sin parámetros compartidos. Los ajustes comparados usan semillas distintas y explícitas.
5. La comparación de zonas se expresa como razón de medianas modeladas; no es una comparación entre todos los segmentos individuales ni una razón de medias aritméticas.
6. Se validan datos positivos y finitos antes de los logaritmos. Las categorías de zona tienen orden explícito y la dimensión de alpha sigue el número de categorías.
7. Se calcula también la probabilidad exacta de superar el umbral urbano usando la cola superior Beta directamente, con lower.tail = FALSE.
8. Metropolis registra explícitamente las aceptaciones y las resume sobre las iteraciones retenidas de las cuatro cadenas. Se rechazan propuestas con log-posterior no finita y se reportan ESS_bulk y ESS_tail.
9. Se agregan diagnósticos modernos, error Monte Carlo, divergencias y alcances de profundidad máxima para los siete ajustes Stan. Los controles de profundidad usan el límite predeterminado de 10 de estas ejecuciones.
10. Se registra sessionInfo() y se actualiza el enlace y contenido en la agenda general.

## Fuentes de diagnóstico

Los criterios y distinciones entre diagnósticos se apoyan en [Stan: diagnóstico de convergencia](https://mc-stan.org/learn-stan/diagnostics-warnings.html) y en la [documentación de posterior](https://mc-stan.org/posterior/reference/ess_bulk.html), consultados al revisar D1S3. Son guías de evaluación del muestreo, no pruebas de validez del modelo.

## Conservación y verificación

El original permanece intacto en Taller/Dia1_S4_Practica. Se mantiene una sola versión de trabajo en esta carpeta. Se reutilizó la caché compilada únicamente cuando el programa Stan era idéntico al ya verificado en D1S2/D1S3; los ajustes y resultados de esta sesión se ejecutan nuevamente.

El detalle de ejecución y compilación queda en resultados_log.txt y render_log.txt. El resultado final de las verificaciones se añade al finalizar.

### Resultado de las verificaciones

- Los siete ajustes Stan y las tres configuraciones del desafío Metropolis se ejecutaron de principio a fin.
- Todos los parámetros revisados en Stan tienen R-hat inferior a 1,004; ESS_bulk superiores a 1400 y ESS_tail superiores a 2000. No se registraron divergencias ni alcances de profundidad máxima.
- Las filas de comparación Bernoulli reportan ahora medias en ambos métodos, con discrepancias compatibles con sus errores Monte Carlo.
- Metropolis ilustra una mala exploración con paso 0,001: aceptación cercana a 0,90, R-hat cercano a 1,148 y ESS_bulk de 18. No se presenta esa configuración como válida para inferencia final. Las otras configuraciones deben evaluarse con sus propios diagnósticos, no solo con su aceptación.
- Avisos de configuración regional C.UTF-8 del entorno no impidieron la ejecución.
- Se separaron láminas de código y de resultados en los ejercicios 1 y 5 para reducir la carga visual. La revisión visual completa de la presentación queda pendiente.

- Compilación final: HTML generado sin avisos de estructura de Quarto. Enlaces del README verificados.


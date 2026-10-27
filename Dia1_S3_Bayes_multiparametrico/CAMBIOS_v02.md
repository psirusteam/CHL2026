# D1S3 · Fundamentos de la inferencia bayesiana II

Versión de trabajo v02 · Revisión del 23 de octubre de 2026.

## Alineación con la nota conceptual

Lunes 23 de noviembre de 2026, 14:30–16:00. Docente: Andrés Gutiérrez. La sesión conserva los cuatro ejes del programa: Normal con media y varianza desconocidas, Multinomial–Dirichlet, simulación MCMC/Metropolis y diagnóstico e introducción a Stan. Incorpora identidad de CEPAL, UNFPA e INE Chile, con logos de CEPAL y UNFPA.

## Revisión conceptual y didáctica

- Conexión explícita con D1S2: dejar de fijar sigma para propagar su incertidumbre.
- Objetivos de aprendizaje y distribución orientativa de los 90 minutos en notas docentes; derivaciones extensas y leapfrog como apoyo para exposición.
- Uso de «parámetros auxiliares» para nuisance, evitando confundirlos con errores del modelo.
- Previa de referencia: impropia, definida respecto de dmu dsigma²; posterior propia para n >= 2 y variación positiva. La t marginal se expresa con su escala al cuadrado explícita.
- Separación entre la solución bajo previa de referencia y el ajuste Stan con previas propias: son modelos distintos, no una verificación exacta de una misma posterior.
- Revisión de las previas amplias en escala logarítmica y la necesidad de evaluar su plausibilidad predictiva.
- Corrección de la independencia posterior: una correlación pequeña entre mu y sigma no implica independencia.
- Multinomial: distinción entre proporciones observadas conocidas y probabilidades bajo el modelo para segmentos comparables.
- Balance detallado, calentamiento y R-hat presentados como partes de la evaluación del muestreo; ninguno garantiza por sí solo una exploración adecuada.
- Fórmula clásica de R-hat identificada como explicación intuitiva; se distingue de la implementación moderna por rangos, división y plegado.
- ESS_bulk y ESS_tail diferenciados; se añaden diagnósticos modernos y controles de divergencias/profundidad del ajuste Normal en Stan.
- Se evitan conclusiones fijas sobre la mejor propuesta Metropolis, una tasa de aceptación universal y la eficiencia por tiempo de cómputo.
- El chequeo predictivo por superposición de densidades se presenta como evaluación parcial, sin declarar validación global del modelo.

## Ajustes del código R

- Comprobaciones de densidad positiva, datos finitos y varianza positiva antes de calcular logaritmos y la solución analítica.
- Cálculo de la previa media Cauchy directamente en logaritmos y rechazo de propuestas con log-posterior no finita.
- Registro explícito de cada aceptación Metropolis; las tasas se calculan con las iteraciones retenidas de las cuatro cadenas, en lugar de inferirlas desde diferencias de una sola cadena.
- Comparación entre algoritmos sobre mu y sigma; se transforma log_sigma antes del ESS de sigma y se busca cada parámetro por nombre en la salida Stan.
- Se conservan semillas explícitas y se registra sessionInfo() junto a los resultados.

## Fuentes para los diagnósticos

- [Manual de Stan: análisis posterior](https://mc-stan.org/docs/2_40/reference-manual/analysis.html).
- [Stan: diagnóstico y resolución de problemas de convergencia](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
- [posterior: ESS_bulk](https://mc-stan.org/posterior/reference/ess_bulk.html).

## Conservación y validación

El original permanece intacto en Taller/Dia1_S3_Bayes_multiparametrico. Esta carpeta contiene una sola versión de trabajo de la sesión. Se actualiza el README general para enlazar esta versión.

La sintaxis R se comprobó antes de ejecutar. El resultado final de ejecución, muestreo y compilación se registra al terminar en esta sección; los detalles quedan en resultados_log.txt y render_log.txt.

### Resultado de la ejecución

- Script completo, dos ajustes Stan y experimentos Metropolis ejecutados correctamente.
- Modelo Normal en Stan: R-hat inferior a 1,004 y ESS_bulk/ESS_tail superiores a 2600 para mu y sigma; sin divergencias ni alcances de profundidad máxima en las cuatro cadenas.
- Metropolis didáctico: mu tiene R-hat 1,001 y ESS_bulk 1598; log_sigma tiene R-hat 1,012 y ESS_bulk 362. Se conserva y se señala explícitamente el diagnóstico insuficiente de log_sigma como parte del ejercicio, sin presentar todas las cadenas como convergidas.
- Aceptación Metropolis: verificada como proporción finita entre cero y uno sobre las iteraciones retenidas de las cuatro cadenas.
- Avisos no bloqueantes: configuración regional C.UTF-8, retirada de metadatos draws_df al resumir funciones y representación de pares con draws combinados. El gráfico de pares ilustra dependencia posterior, no diagnóstico entre cadenas.
- Se inspeccionó el gráfico actualizado de posterior conjunta y marginales. La revisión visual completa de las diapositivas queda pendiente.

- Compilación final con Quarto: HTML generado correctamente con todos los resultados actualizados.


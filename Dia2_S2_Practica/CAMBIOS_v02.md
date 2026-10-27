# D2S2 · Práctica computacional: modelos lineales bayesianos

Versión de trabajo v02 · Revisión del 23 de octubre de 2026.

## Programa y organización

Martes 24 de noviembre, 11:30–13:00. Docente: Stalyn Guerrero. Continúa D2S1 con comparación de regiones, evaluación predictiva, dispersión por zona y chequeo de previas. El ejercicio de población implícita se conserva como extensión y puente al modelo de conteos. Se propone un recorrido de 90 minutos y un producto por pareja.

## Cambios principales

- Identidad institucional, fecha, horario y logos coherentes con la nota conceptual.
- Ajustes regionales conservados como objetos Stan para diagnosticarlos; columnas de draws nombradas explícitamente en lugar de depender de un renombrado por posición de toda la salida.
- Semillas distintas por ajuste y calentamiento explícito. Se reutiliza el ajuste completo de Aysén como modelo B porque datos, diseño y previas son idénticos, evitando ajustar dos veces el mismo modelo.
- Distribución posterior de la diferencia Magallanes menos Aysén, con intervalos y probabilidad positiva. Un cambio de signo de medias no demuestra diferencia ni determina una estructura nacional.
- PSIS-LOO incorpora relative_eff sobre la verosimilitud y los identificadores de cadena. La verosimilitud se reescala por observación para evitar problemas de exponenciación; ese reescalado no cambia su eficiencia relativa.
- Diferencia elpd B menos A calculada por nombre y contribuciones puntuales, sin asumir el orden que devuelve loo_compare. Su error estándar se calcula con las diferencias emparejadas por observación.
- Diagnósticos Pareto-k para A, B y C con umbral dependiente del número de draws. Los rankings se presentan como provisionales cuando hay observaciones señaladas, hasta comprobar la aproximación con reajustes exactos o una estrategia alternativa de validación.
- Se distingue LOO entre segmentos aptos de una región de la predicción en otras regiones, no aptos o bloques espaciales. p_loo no se equipara automáticamente al número de coeficientes.
- Dispersión por zona evaluada con un chequeo parcial, sin declarar que un único estadístico valida todo el modelo.
- Los totales del ejercicio 4 son cantidades continuas implícitas condicionadas a viviendas observadas; no incorporan ruido Poisson de conteos ni constituyen una estimación validada de omisión censal. La corrección V/c requiere cobertura válida positiva; con cobertura cero se necesita otra fuente de exposición.
- Chequeo previo calculado en log-densidad para conservar extremos sin desbordamiento. Se reportan porcentajes que desbordarían al exponentiar y que quedan fuera del gráfico. Las líneas 1 y 6 se identifican como referencias ilustrativas, no límites universales.
- Notación de las previas alternativas con varianzas explícitas, coherente con las desviaciones estándar del código.
- Diagnósticos modernos de ambos ajustes regionales, modelo reducido y modelo con sigma por zona; controles de NUTS y sessionInfo().

## Referencias técnicas consultadas

- [loo: relative_eff](https://mc-stan.org/loo/reference/relative_eff.html).
- [Programas Stan para loo](https://mc-stan.org/loo/articles/loo2-with-rstan.html).
- [loo: estimación y comparación predictiva](https://mc-stan.org/loo/reference/loo.html).

## Conservación y verificación

El material original permanece intacto en Taller/Dia2_S2_Practica. Se mantiene una sola versión de trabajo. Se reutilizó la caché compilada del programa lineal únicamente tras comprobar que su fuente Stan era idéntica a la de D2S1; los ajustes de esta sesión se ejecutaron nuevamente.

Los detalles se registran en resultados_log.txt y render_log.txt. Las verificaciones finales se añaden al terminar.

## Uso y codificación

El código y los documentos usan UTF-8. En RStudio se recomienda abrirlos con esa codificación. Durante la verificación en este entorno Windows, la configuración inicial C.UTF-8 no era válida; se estableció LC_CTYPE en English_United States.utf8 antes de ejecutar source(..., encoding = "UTF-8") y regenerar los gráficos. Este ajuste corrigió acentos y símbolos sin modificar la metodología.

## Límite de la comparación predictiva

En la ejecución verificada, PSIS señaló una observación en B (Pareto-k superior a 1) y una en C (superior a 0,7). Se identifican expresamente en la presentación. Los valores de elpd se conservan como salida del ejercicio, pero su ranking no se adopta como conclusión validada: requiere reajuste exacto de los casos señalados o una validación alternativa. Reconocer y tratar este diagnóstico forma parte de la práctica.

### Resultado final de las verificaciones

- Ejecución completa de los cuatro ajustes y de todos los ejercicios finalizada correctamente, con regeneración final en UTF-8.
- R-hat máximo 1,0063; ESS_bulk mínimo superior a 1500 y ESS_tail mínimo superior a 1800. Sin divergencias ni alcances de profundidad máxima en los ajustes revisados.
- PSIS-LOO: A sin casos señalados; B con un caso de Pareto-k 1,061; C con un caso de Pareto-k 0,838. Las conclusiones predictivas quedan condicionadas a comprobar estas aproximaciones.
- El chequeo previo conserva y cuantifica los extremos en escala logarítmica. Se verificó que los porcentajes de desbordamiento fueran finitos y estuvieran entre cero y cien.
- Se inspeccionó el gráfico del chequeo previo regenerado y se confirmó la corrección de acentos y símbolos. La revisión visual completa de las diapositivas queda pendiente.
- HTML compilado con todos los resultados actualizados. Los avisos PSIS se conservan como diagnósticos del ejercicio, sin ocultarlos ni presentar el ranking como validado.

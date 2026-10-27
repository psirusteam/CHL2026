# D1S2 · Fundamentos de la inferencia bayesiana I

Versión de trabajo v02, 21 de octubre de 2026.

## Correspondencia con la nota conceptual

Sesión del lunes 23 de noviembre, 11:30–13:00. Docente: Andrés Gutiérrez. Introduce los modelos uniparamétricos como base para la posterior especificación del modelo de población. Conserva los tres ejemplos del material inicial: Bernoulli–Beta, Poisson–Gamma con exposición y Normal con varianza fijada.

## Ajustes realizados

- Título institucional, docente, fecha y horario alineados con la nueva nota conceptual.
- Logos de CEPAL y UNFPA incorporados en la presentación HTML.
- Objetivos de aprendizaje y distribución orientativa de los 90 minutos en notas docentes. La ejecución guiada de Stan se retoma en la práctica de la tarde.
- Distinción entre proporción descriptiva conocida del conjunto de segmentos y parámetro Bernoulli para segmentos comparables.
- Corrección de la asociación entre pérdida 0–1 y MAP para parámetros continuos.
- Precisión de la definición de conjunto de máxima densidad, que puede ser disjunto.
- Revisión de las interpretaciones de las previas Beta y Gamma; explicitación de la parametrización shape–rate.
- Explicación del uso didáctico de la desviación observada como sigma fija y de la incertidumbre que se omite.
- Sustitución de conclusiones numéricas fijas por interpretaciones dependientes de las tablas y los resultados.
- Cálculo de la Beta-Binomial en escala logarítmica para evitar cocientes de funciones beta que pueden desbordarse o reducirse numéricamente a cero.
- Corrección de las etiquetas del gráfico de contracción: exp(media de log-densidad) no es la media aritmética de personas por vivienda.

## Conservación de materiales

La carpeta inicial Taller/Dia1_S2_Bayes_uniparametrico permanece intacta. Aquí se mantiene una sola versión de trabajo, junto con sus datos, código, modelos Stan, resultados y recursos gráficos. No se generan respaldos intermedios dentro de esta carpeta.

## Validación

Se revisaron los tres programas Stan y la correspondencia de las fórmulas conjugadas con el código. La ejecución y compilación se registran en resultados_log.txt y render_log.txt, respectivamente. El estado final de estas verificaciones se actualiza al terminar.

## Aspectos didácticos para la revisión conjunta

La cantidad de derivaciones y ejemplos computacionales es exigente para 90 minutos. Se conserva todo el contenido disponible: las notas docentes proponen priorizar la interpretación, los dos primeros modelos y la conexión con la exposición, usando las derivaciones extensas como apoyo. La reubicación definitiva de material se decidirá al revisar la sesión.

### Resultado de las verificaciones

- Ejecución completa del script R y de los tres ajustes Stan: finalizada correctamente.
- Predictiva Beta-Binomial: probabilidades finitas y suma igual a uno con tolerancia de 1e-8.
- R-hat de los parámetros theta, delta y mu: entre 1,0007 y 1,0038; tamaños efectivos reportados por RStan superiores a 1100. Estos indicadores no sustituyen la revisión de trazas y del ajuste del modelo.
- Compilación final con Quarto 1.9.36: HTML generado con los resultados actualizados.
- Se observaron avisos de configuración regional C.UTF-8 del entorno; no impidieron ejecutar los análisis ni compilar.
- Pendiente revisión visual conjunta de la presentación y decisión sobre la profundidad de las derivaciones para los 90 minutos disponibles.

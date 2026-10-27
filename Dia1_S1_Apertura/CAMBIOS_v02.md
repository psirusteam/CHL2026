# Revisión de Día 1 · Sesión 1: Apertura

Versión: v02 · 21 de octubre de 2026 · Propuesta para revisión conjunta.
Fuente rectora: nota_conceptual_taller_v02_revAG_revSJ_revAG2.docx.
Original: Taller/Dia1_S1_Apertura, conservado sin modificaciones.

## Cambios realizados

| Elemento | Material original | Versión 02 |
|---|---|---|
| Título | Modelación bayesiana para estimación de población | Título de la nueva nota conceptual |
| Instituciones | Predominio de CEPAL | CEPAL, UNFPA e INE Chile |
| Propósito | Actualización de proyecciones con «horizonte 2026» | Evaluación complementaria de omisión y reconciliación para proyecciones nacionales y subnacionales |
| Objetivos | Incluía un sexto objetivo adicional | Cinco objetivos específicos y objetivo general de la nota |
| Fechas y horario | Agenda desde 09:00; segundo bloque 11:00–13:00 | 23–26 noviembre; 09:30–11:00 y 11:30–13:00 |
| Miércoles | Conteos, variables, implementación | Variables satelitales I y II; modelo de Chile; práctica |
| Jueves | Modelo de Chile; diagnóstico y predicción | Diagnóstico; predicción; práctica; cierre institucional |
| Apertura y cierre | Sin intervenciones institucionales detalladas | Participantes de la nota conceptual |
| Profundidad | Desarrollo matemático extenso en apertura | Introducción conceptual; desarrollo anterior conservado como complemento |

## Conservación de información

- La carpeta original permanece intacta y constituye la versión inicial (v01).
- El QMD original en Taller/Dia1_S1_Apertura conserva cobertura, tablas, gráficos, Bayes, fórmulas del modelo, previas, predicción, agregación y notación. No se duplica dentro de la versión de trabajo.
- Se copiaron Rcode, data, resultados, img, referencias y estilos, sin cambios de contenido.
- No se copió el HTML anterior como salida de v02, para evitar confundirlo con la presentación revisada.
- Las librerías generadas de la presentación anterior no son fuentes; hubo errores de copia en esas librerías por rutas extensas en la unidad de red. No se necesitan para generar nuevamente la presentación.

## Temas que requieren revisión en las sesiones posteriores

1. Reubicar el contenido de conteos del antiguo miércoles por la mañana en el martes y/o la sesión del modelo de Chile, según su extensión.
2. Desarrollar las dos sesiones de variables satelitales con los materiales del docente indicado en la nota, sin suponer que el material actual cubre ambas.
3. Integrar la implementación y el modelo de Chile en el miércoles por la tarde.
4. Separar diagnóstico y predicción en las dos sesiones del jueves.
5. Contrastar con los datos y el código la definición de cobertura, el tratamiento de cobertura cero, la exposición corregida, el techo de densidad, las previas y el efecto regional. La nota conceptual no valida esos detalles operativos.
6. Revisar las afirmaciones generales sobre la ronda censal y la comparación frecuentista/bayesiana antes de reutilizar el complemento en una sesión técnica.

## Validación

Se verificó el contenido de la apertura revisada contra el texto y el cronograma de la nota conceptual. La apertura nueva es estática: no modifica análisis ni necesita ejecutar R. Queda pendiente compilar y revisar visualmente la presentación con Quarto, que no se encontró en el entorno disponible. No se presenta esta versión como HTML final validado.

## Uso

Desde esta carpeta, con Quarto instalado:

    quarto render Apertura.qmd

La presentación principal tiene 17 láminas de contenido, además de la portada. Los tiempos internos de bienvenida y presentación técnica no están definidos por la nota; se dejan a la coordinación docente.


## Incorporación de identidad UNFPA

Se agregó el logo vectorial existente en doc/UNFPA_logo.svg junto al de CEPAL, mediante logos_institucionales.html y estilos compartidos para portada y diapositivas en HTML. Durante los ajustes se mantiene una sola versión de trabajo; se eliminaron los respaldos intermedios redundantes. Pendiente compilación y revisión visual con Quarto. Esta franja está destinada a la presentación HTML y se oculta en impresión.

## Criterio de trabajo acordado

Mantener una sola versión de trabajo por sesión durante los ajustes. Actualizar esta misma versión y su registro de cambios, sin generar copias intermedias. Conservar el material original en su ubicación inicial.


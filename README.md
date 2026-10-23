# Taller de fortalecimiento de capacidades en modelación estadística de población en Chile

**CEPAL · Fondo de Población de las Naciones Unidas (UNFPA) · Instituto Nacional de Estadísticas de Chile (INE)**

**Fecha:** 23 al 26 de noviembre de 2026.  
**Modalidad:** presencial.  
**Lugar:** sala Carmen Miró, CEPAL, Santiago de Chile.  
**Horario:** 09:30 a 18:00.

## Introducción

Durante 2026, la CEPAL y UNFPA-LACRO han acompañado al INE de Chile en el desarrollo de un modelo de población por segmento censal orientado a la estimación de la omisión censal. La metodología combina información censal y covariables satelitales mediante un modelo jerárquico Poisson–lognormal ajustado con inferencia bayesiana.

El taller busca fortalecer las capacidades del equipo técnico del INE para especificar, implementar, diagnosticar y validar estos modelos, como herramienta complementaria para la evaluación de la omisión censal y la reconciliación de las estimaciones de población que sirven de base para la actualización de las proyecciones nacionales y subnacionales.

Las sesiones combinan exposiciones teóricas con prácticas computacionales guiadas en R y Stan, utilizando datos del caso chileno. El recorrido parte de los fundamentos de la inferencia bayesiana y avanza hacia la incorporación de covariables satelitales, la implementación del modelo de población, su diagnóstico y la generación de predicciones por segmento censal.

Al finalizar, se espera que el equipo pueda reproducir el flujo completo de trabajo y disponga de códigos, bases y documentación para continuar la aplicación y adaptación institucional de la metodología.

## Agenda

| Día | Fecha | Sesión | Contenido |
|---|---|---|---|
| Día 1 | 23/11/2026 | [Apertura](Dia1_S1_Apertura/Apertura.qmd) | 1. Bienvenida de CEPAL, UNFPA e INE Chile<br>2. Contexto de la cooperación técnica y propósito del modelo<br>3. Objetivos, metodología y resultados esperados<br>4. Organización y materiales del taller |
| Día 1 | 23/11/2026 | [Fundamentos de la inferencia bayesiana I: modelos uniparamétricos](../../Dia1_S2_Bayes_uniparametrico/Bayes_uniparametrico.qmd) | 1. Regla de Bayes, distribución previa y posterior<br>2. Modelo Bernoulli–Beta<br>3. Modelo Poisson–Gamma con exposición<br>4. Modelo Normal y chequeo predictivo |
| Día 1 | 23/11/2026 | [Fundamentos de la inferencia bayesiana II: modelos multiparamétricos](../../Dia1_S3_Bayes_multiparametrico/Bayes_multiparametrico.qmd) | 1. Modelo Normal con media y varianza desconocidas<br>2. Modelo Multinomial–Dirichlet<br>3. Simulación MCMC y algoritmo de Metropolis<br>4. Diagnósticos de convergencia e introducción a Stan |
| Día 1 | 23/11/2026 | [Práctica computacional en RStan](../../Dia1_S4_Practica/Practica_Dia1.qmd) | 1. Implementación de modelos uniparamétricos y multiparamétricos<br>2. Ejercicios con datos de Magallanes<br>3. Comparación entre regiones y zonas<br>4. Exploración del algoritmo de Metropolis |
| Día 2 | 24/11/2026 | [Introducción a los modelos lineales con métodos bayesianos](../../Dia2_S1_Modelo_lineal/Modelo_lineal.qmd) | 1. Regresión bayesiana del logaritmo de personas por vivienda<br>2. Incorporación e interpretación de covariables<br>3. Comparación con mínimos cuadrados ordinarios<br>4. Ajuste, chequeo predictivo y limitaciones |
| Día 2 | 24/11/2026 | [Práctica computacional en RStan](../../Dia2_S2_Practica/Practica_Dia2_lineal.qmd) | 1. Ajuste de modelos lineales<br>2. Efectos por región y dispersión por zona<br>3. Comparación mediante LOO<br>4. Predicción y chequeo predictivo previo |
| Día 2 | 24/11/2026 | [Modelo Poisson con vínculo loglineal](../../Dia2_S3_Modelo_Poisson/Modelo_Poisson.qmd) | 1. Modelo Poisson para conteos<br>2. Offset de viviendas particulares ocupadas<br>3. Sobredispersión y modelo Poisson–lognormal<br>4. Estructura jerárquica y encogimiento |
| Día 2 | 24/11/2026 | [Práctica computacional en RStan](../../Dia2_S4_Practica/Practica_Dia2_Poisson.qmd) | 1. Ajuste y comparación de modelos de conteos<br>2. Modelo binomial negativo<br>3. Predicción de totales regionales<br>4. Parametrizaciones centrada y no centrada |
| Día 3 | 25/11/2026 | Variables satelitales para modelos de población I | 1. Fuentes de información satelital y geoespacial<br>2. Cobertura y resolución de los insumos<br>3. Vinculación con segmentos censales<br>4. Estadísticas zonales |
| Día 3 | 25/11/2026 | Variables satelitales para modelos de población II | 1. Preparación, transformación y estandarización de covariables<br>2. Relación con las personas por vivienda<br>3. Colinealidad y selección de predictores<br>4. Integración de covariables en el modelo de población |
| Día 3 | 25/11/2026 | El modelo de población en Chile | 1. Especificación jerárquica Poisson–lognormal<br>2. Elegibilidad, cobertura y exposición<br>3. Datos e implementación en Stan<br>4. Flujo de ajuste del modelo por región |
| Día 3 | 25/11/2026 | [Práctica computacional en RStan](../../Dia3_S4_Practica/Practica_Dia3.qmd) | 1. Ajuste del modelo en Magallanes<br>2. Generación de predicciones<br>3. Escenarios de cobertura<br>4. Exploración de especificaciones alternativas |
| Día 4 | 26/11/2026 | Diagnóstico del modelo | 1. Convergencia: R-hat y tamaño efectivo de muestra<br>2. Trazas, rangos y relaciones entre parámetros<br>3. Identificabilidad y problemas de muestreo<br>4. Chequeo predictivo y validación del ajuste |
| Día 4 | 26/11/2026 | Predicción en los segmentos censales | 1. Predicción a partir de las cadenas posteriores<br>2. Incorporación de escenarios de cobertura<br>3. Agregación por comuna y región<br>4. Cuantificación de incertidumbre y contraste de resultados |
| Día 4 | 26/11/2026 | [Práctica computacional en RStan](../../Dia4_S3_Practica/Practica_Dia4.qmd) | 1. Diagnóstico de las cadenas del modelo de Chile<br>2. Revisión de la especificación<br>3. Construcción de escenarios de cobertura<br>4. Interpretación de resultados y discusión de validación |
| Día 4 | 26/11/2026 | Cierre del evento | 1. Balance del taller<br>2. Aplicación y adaptación institucional de la metodología<br>3. Cierre de CEPAL y UNFPA |

Las fechas y sesiones siguen la nueva nota conceptual. Los contenidos constituyen una propuesta de organización a partir de los materiales disponibles y se ajustarán al revisar cada sesión. Los enlaces apuntan a la versión de trabajo de Apertura y a los materiales originales de las sesiones aún no revisadas. Las sesiones reorganizadas se enlazarán cuando se prepare su versión correspondiente.

**Horario diario:** primer bloque de 09:30 a 11:00; segundo de 11:30 a 13:00; tercero de 14:30 a 16:00; cuarto de 16:30 a 18:00. Pausas de café de 11:00 a 11:30 y de 16:00 a 16:30; almuerzo de 13:00 a 14:30.

## Participantes y requisitos

El taller está dirigido al personal técnico del INE de Chile vinculado a la estimación y proyección de población, en particular a los equipos responsables del procesamiento de información censal.

Se recomiendan conocimientos sólidos de estadística —manejo de bases de datos, inferencia y regresión— y manejo de R y RStudio. No se requiere experiencia previa en Stan ni en modelación bayesiana. Se solicita contar con computador con conexión a internet y con R, RStudio y Stan previamente instalados.

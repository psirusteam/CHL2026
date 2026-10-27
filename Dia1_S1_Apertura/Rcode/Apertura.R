##############################################################################
## Dia 1 - Sesion 1: Apertura                                               ##
## Resultados que usa la presentacion Apertura.qmd.                         ##
## Ejecutar con el directorio de trabajo en la carpeta Dia1_S1_Apertura/:   ##
##   source("Rcode/Apertura.R")                                             ##
##############################################################################

source("Rcode/utilidades.R")

## ---- paquetes ----
library(tidyverse)
theme_set(theme_minimal(base_size = 15))

## ---- datos ----
resumen <- readRDS("data/resumen_nacional_segmentos.rds")
aysen <- readRDS("data/segmentos_aysen.rds")
magallanes <- readRDS("data/segmentos_magallanes.rds")
diccionario <- readRDS("data/diccionario_covariables.rds")

## ---- tabla-nacional ----
tabla_nacional <- resumen %>%
  group_by(clasificacion = clasificacion_cobertura) %>%
  summarise(segmentos = sum(n_segmentos),
            aptos = sum(n_segmentos[apto_modelo]),
            personas = sum(n_per),
            viviendas_ocupadas = sum(n_vpo)) %>%
  mutate(pct_segmentos = 100 * segmentos / sum(segmentos))
tabla_nacional

## ---- guardar-tabla-nacional ----
guardar(tabla_nacional, "tabla_nacional")
guardar_valor("segmentos_pais", sum(tabla_nacional$segmentos))
guardar_valor("personas_pais", sum(tabla_nacional$personas))

## ---- cobertura-region ----
grafico_region <- resumen %>%
  group_by(cod_region, clasificacion_cobertura) %>%
  summarise(n = sum(n_segmentos), .groups = "drop") %>%
  group_by(cod_region) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(factor(cod_region), pct, fill = clasificacion_cobertura)) +
  geom_col() +
  scale_fill_manual(values = c(Completa = "#2e75b6", Parcial = "#f4b183", Incompleta = "#c00000")) +
  labs(x = "Región (código)", y = "% de segmentos", fill = "Cobertura")
guardar_figura(grafico_region, "grafico_region", alto = 5)

## ---- densidad-clasificacion ----
tabla_densidad <- resumen %>%
  filter(apto_modelo) %>%
  group_by(zona, clasificacion_cobertura) %>%
  summarise(delta = sum(n_per) / sum(n_vpo), .groups = "drop") %>%
  pivot_wider(names_from = clasificacion_cobertura, values_from = delta)
guardar(tabla_densidad, "tabla_densidad")

## ---- ejemplo-bayes ----
# Densidad de Aysen (personas por vivienda) en segmentos de cobertura completa:
# previa Gamma(a, b) y posterior Gamma(a + sum(Y), b + sum(V))
completos <- aysen %>% filter(clasificacion_cobertura == "Completa", !excluido_colectivo, n_vpo > 0)
a <- 30; b <- 10                                   # previa: media 3, desviacion ~0,55
a_post <- a + sum(completos$n_per)
b_post <- b + sum(completos$n_vpo)

ejemplo <- tibble(
  momento = c("Previa", "Posterior"),
  media = c(a / b, a_post / b_post),
  desv_est = c(sqrt(a) / b, sqrt(a_post) / b_post),
  q2.5 = c(qgamma(0.025, a, b), qgamma(0.025, a_post, b_post)),
  q97.5 = c(qgamma(0.975, a, b), qgamma(0.975, a_post, b_post))
)
ejemplo

## ---- guardar-ejemplo-bayes ----
guardar(ejemplo, "ejemplo_bayes")
guardar_valor("n_completos_aysen", nrow(completos))
guardar_valor("personas_completos_aysen", sum(completos$n_per))
guardar_valor("viviendas_completos_aysen", sum(completos$n_vpo))

grafico_ejemplo <- tibble(delta = seq(1.5, 4.5, length.out = 500)) %>%
  mutate(Previa = dgamma(delta, a, b), Posterior = dgamma(delta, a_post, b_post)) %>%
  pivot_longer(-delta, names_to = "distribucion", values_to = "densidad") %>%
  mutate(distribucion = factor(distribucion, levels = c("Previa", "Posterior"))) %>%
  ggplot(aes(delta, densidad, colour = distribucion)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ distribucion, scales = "free_y") +
  scale_colour_manual(values = c(Previa = "#7f7f7f", Posterior = "#1f4e79"), guide = "none") +
  labs(x = "δ: personas por vivienda ocupada", y = "Densidad")
guardar_figura(grafico_ejemplo, "grafico_ejemplo", alto = 4.2)

## ---- datos-taller ----
tabla_regiones <- bind_rows(Aysén = aysen, Magallanes = magallanes, .id = "region") %>%
  group_by(region) %>%
  summarise(segmentos = n(), aptos = sum(apto_modelo),
            completos = sum(clasificacion_cobertura == "Completa"),
            con_colectivas = sum(excluido_colectivo),
            personas = sum(n_per), viviendas_ocupadas = sum(n_vpo),
            personas_por_vivienda = sum(n_per) / sum(n_vpo))
tabla_regiones

## ---- guardar-datos-taller ----
guardar(tabla_regiones, "tabla_regiones")

## ---- covariables-taller ----
covariables <- sub("_std$", "", grep("_std$", names(aysen), value = TRUE))
tabla_covariables <- diccionario %>%
  filter(variable %in% covariables) %>%
  select(variable, fuente, descripcion = descripcion_es) %>%
  mutate(descripcion = str_trunc(descripcion, 80))
guardar(tabla_covariables, "tabla_covariables")

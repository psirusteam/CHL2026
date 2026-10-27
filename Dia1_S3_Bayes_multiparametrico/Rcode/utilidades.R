##############################################################################
## Utilidades comunes del taller                                            ##
##                                                                          ##
## Los scripts de Rcode/ corren los analisis y GUARDAN sus resultados en    ##
## resultados/ (tablas y textos .rds, figuras .png, valores numericos).     ##
## Las presentaciones (.qmd) no ejecutan codigo pesado: muestran el codigo  ##
## de los scripts y LEEN los resultados guardados con mostrar().            ##
##                                                                          ##
## Todas las rutas son relativas a la carpeta de la sesion.                 ##
##############################################################################

carpeta_resultados <- "resultados"

ruta_resultado <- function(nombre, ext = "rds") {
  file.path(carpeta_resultados, paste0(nombre, ".", ext))
}

# Guardar una tabla (data.frame), un vector o cualquier objeto de R
guardar <- function(objeto, nombre) {
  dir.create(carpeta_resultados, showWarnings = FALSE, recursive = TRUE)
  saveRDS(objeto, ruta_resultado(nombre))
  invisible(objeto)
}

# Guardar la salida impresa de un objeto (lo que se veria en la consola)
guardar_texto <- function(objeto, nombre) {
  guardar(structure(capture.output(print(objeto)), class = "salida_consola"), nombre)
  invisible(objeto)
}

# Guardar un grafico de ggplot (o bayesplot / patchwork) como PNG
guardar_figura <- function(grafico, nombre, ancho = 12, alto = 5.5) {
  dir.create(carpeta_resultados, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(ruta_resultado(nombre, "png"), grafico, width = ancho, height = alto,
                  dpi = 110, bg = "white")
  invisible(grafico)
}

# Valores numericos que se citan en el texto de la presentacion
valores <- list()
guardar_valor <- function(nombre, valor) {
  valores[[nombre]] <<- valor
  guardar(valores, "valores")
  invisible(valor)
}

# Tabla con formato de la CEPAL (separador de miles "." y decimal ",")
tba <- function(dat, cap = NA, digits = 3, font_size = NULL) {
  if (is.data.frame(dat)) {
    dat <- dplyr::mutate(dat, dplyr::across(dplyr::any_of(c("ess_bulk", "ess_tail", "n_eff")), round))
  }
  tabla <- knitr::kable(dat, format = "html", digits = digits, caption = cap,
                        format.args = list(big.mark = ".", decimal.mark = ",")) |>
    kableExtra::kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE,
                              font_size = font_size) |>
    kableExtra::kable_classic(full_width = FALSE, html_font = "Arial Narrow")
  tabla
}

# Mostrar en la presentacion un resultado guardado: figura, tabla o texto
mostrar <- function(nombre, ...) {
  png <- ruta_resultado(nombre, "png")
  if (file.exists(png)) return(knitr::include_graphics(png))
  archivo <- ruta_resultado(nombre)
  if (!file.exists(archivo)) {
    return(cat("[Resultado '", nombre, "' no encontrado: ejecute primero el script de Rcode/]\n", sep = ""))
  }
  x <- readRDS(archivo)
  if (inherits(x, "salida_consola")) {
    cat(unclass(x), sep = "\n")
    return(invisible(NULL))
  }
  if (is.data.frame(x) || is.matrix(x)) return(tba(x, ...))
  print(x)
  invisible(NULL)
}

# Valores guardados para citar en el texto
leer_valores <- function() {
  if (file.exists(ruta_resultado("valores"))) readRDS(ruta_resultado("valores")) else list()
}

# Formato numerico para el texto
fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")

# Resumen posterior de un ajuste de rstan como data.frame
resumen_stan <- function(fit, pars) {
  s <- rstan::summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  data.frame(parametro = rownames(s), media = s[, "mean"], desv_est = s[, "sd"],
             q2.5 = s[, "2.5%"], mediana = s[, "50%"], q97.5 = s[, "97.5%"],
             n_eff = round(s[, "n_eff"]), Rhat = s[, "Rhat"], row.names = NULL)
}

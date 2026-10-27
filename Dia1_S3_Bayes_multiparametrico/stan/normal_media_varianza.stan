// Modelo Normal con media y varianza desconocidas: log(personas por vivienda)
data {
  int<lower=0> n;
  vector[n] y;                 // log(delta_i)
}
parameters {
  real mu;                     // media de log(delta)
  real<lower=0> sigma;         // desviacion estandar de log(delta)
}
model {
  mu ~ normal(0, 10);          // previas debilmente informativas
  sigma ~ cauchy(0, 2.5);      // media Cauchy positiva por la restriccion lower=0
  y ~ normal(mu, sigma);
}
generated quantities {
  // Personas por vivienda en un segmento nuevo
  real delta_nuevo = exp(normal_rng(mu, sigma));
}

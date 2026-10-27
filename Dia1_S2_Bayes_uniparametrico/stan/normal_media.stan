// Modelo Normal con varianza conocida: media de log(delta)
data {
  int<lower=0> n;
  vector[n] y;                 // log(personas por vivienda)
  real<lower=0> sigma;         // desviacion estandar conocida
  real mu0;                    // previa N(mu0, tau0^2)
  real<lower=0> tau0;
}
parameters {
  real mu;
}
model {
  mu ~ normal(mu0, tau0);
  y ~ normal(mu, sigma);
}

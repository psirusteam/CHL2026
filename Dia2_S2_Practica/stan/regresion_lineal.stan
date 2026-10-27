// Regresion lineal bayesiana para log(personas por vivienda)
data {
  int<lower=1> n;              // segmentos
  int<lower=1> p;              // columnas de X (incluye el intercepto)
  matrix[n, p] X;              // covariables
  vector[n] y;                 // log(delta_i)
  real<lower=0> sd_beta;       // desviacion estandar de la previa de beta
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;
}
model {
  beta ~ normal(0, sd_beta);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(X * beta, sigma);
}
generated quantities {
  // Replicas de los datos para el chequeo predictivo
  array[n] real y_rep = normal_rng(X * beta, sigma);
}

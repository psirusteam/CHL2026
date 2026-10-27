// Regresion lineal con una desviacion estandar distinta por zona
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  vector[n] y;                               // log(delta_i)
  int<lower=1> J;                            // numero de zonas
  array[n] int<lower=1, upper=J> zona;       // zona de cada segmento
}
parameters {
  vector[p] beta;
  vector<lower=0>[J] sigma;                  // una sigma por zona
}
model {
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(X * beta, sigma[zona]);         // indexacion: sigma de la zona de cada segmento
}
generated quantities {
  array[n] real y_rep = normal_rng(X * beta, sigma[zona]);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | X[i] * beta, sigma[zona[i]]);
}

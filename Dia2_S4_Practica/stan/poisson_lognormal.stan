// Poisson-lognormal por zona, con la misma exposicion del Poisson de la practica.
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  array[n] int<lower=0> y;
  vector[n] log_v;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;
  vector[n] z;
}
model {
  beta ~ normal(0, 1);
  sigma ~ normal(0, 0.5);
  z ~ std_normal();
  y ~ poisson_log(log_v + X * beta + sigma * z);
}
generated quantities {
  array[n] int y_rep;
  int<lower=0> rng_out_of_range = 0;
  // Cada replica genera un efecto nuevo para el mismo diseno observado.
  for (i in 1:n) {
    real log_lambda = log_v[i] + X[i] * beta + sigma * normal_rng(0, 1);
    if (log_lambda >= 30 * log(2)) {
      y_rep[i] = -1;
      rng_out_of_range += 1;
    } else {
      y_rep[i] = poisson_log_rng(log_lambda);
    }
  }
}

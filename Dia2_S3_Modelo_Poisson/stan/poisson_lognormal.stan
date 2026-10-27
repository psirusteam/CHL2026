// Poisson-lognormal: un efecto aleatorio por segmento captura la sobredispersion
//   Y_i ~ Poisson(V_i * delta_i),  log(delta_i) = x_i' beta + sigma * z_i,  z_i ~ N(0, 1)
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  array[n] int<lower=0> y;
  vector[n] log_v;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;         // variabilidad de log(delta) entre segmentos
  vector[n] z;                 // efecto estandarizado de cada segmento
}
transformed parameters {
  vector[n] log_delta = X * beta + sigma * z;
}
model {
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  z ~ std_normal();
  y ~ poisson_log(log_v + log_delta);
}
generated quantities {
  // Mixed replication: new latent density at each observed design/exposure.
  array[n] int y_rep;
  int<lower=0> rng_out_of_range = 0;
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

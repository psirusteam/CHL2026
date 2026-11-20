// Poisson-lognormal: un efecto aleatorio por segmento captura la sobredispersion
//   Y_i ~ Poisson(mu_i),  log(mu_i) = x_i' beta + sigma * z_i,  z_i ~ N(0, 1)
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  array[n] int<lower=0> y;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma;         // variabilidad de log(mu) entre segmentos
  vector[n] z;                 // efecto estandarizado de cada segmento
}
transformed parameters {
  vector[n] log_mu = X * beta + sigma * z;
}
model {
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  z ~ std_normal();
  y ~ poisson_log(log_mu);
}
generated quantities {
  // Mixed replication: new latent effect at each observed design.
  array[n] int y_rep;
  int<lower=0> rng_out_of_range = 0;
  for (i in 1:n) {
    real log_lambda = X[i] * beta + sigma * normal_rng(0, 1);
    if (log_lambda >= 30 * log(2)) {
      y_rep[i] = -1;
      rng_out_of_range += 1;
    } else {
      y_rep[i] = poisson_log_rng(log_lambda);
    }
  }
}

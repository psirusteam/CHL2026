// Regresion Poisson para conteos: log(lambda) = X * beta
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  array[n] int<lower=0> y;     // personas
}
parameters {
  vector[p] beta;
}
model {
  beta ~ normal(0, 1);
  y ~ poisson_log(X * beta);
}
generated quantities {
  array[n] int y_rep;
  int<lower=0> rng_out_of_range = 0;
  vector[n] log_lik;
  for (i in 1:n) {
    real log_lambda = X[i] * beta;
    if (log_lambda >= 30 * log(2)) {
      y_rep[i] = -1; // invalid simulation marker, never a population count
      rng_out_of_range += 1;
    } else {
      y_rep[i] = poisson_log_rng(log_lambda);
    }
    log_lik[i] = poisson_log_lpmf(y[i] | log_lambda);
  }
}

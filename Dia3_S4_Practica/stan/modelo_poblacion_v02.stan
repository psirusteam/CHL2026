// Modelo Poisson-lognormal para conteos poblacionales por segmento censal.
// Version 02: incluye la prediccion completa (D_tot: X_tot, Z_tot, V_tot,
// Y_tot, densidad_hat) dentro del propio programa Stan. 
//

data {
  int<lower=1> D_obs;                  // Segmentos aptos para el modelo (S_obs), de la region
  int<lower=1> D_tot;                  // Total de segmentos de la region
  int<lower=1> P;                      // Columnas de X (efectos fijos)
  int<lower=1> Q;                      // Columnas de Z (1 en el ajuste por region)
  array[D_obs] int<lower=0> Y_obs;     // Personas censadas (segmentos aptos)
  array[D_obs] int<lower=0> V_obs;     // Viviendas particulares ocupadas (offset, aptos)
  array[D_tot] int<lower=0> V_tot;     // Offset para todos los segmentos (censo o satelital)
  matrix[D_obs, P] X_obs;              // Covariables, segmentos aptos
  matrix[D_tot, P] X_tot;              // Covariables, todos los segmentos
  matrix[D_obs, Q] Z_obs;              // Indicadores de region, segmentos aptos
  matrix[D_tot, Q] Z_tot;              // Indicadores de region, todos los segmentos
  real<lower=0> delta_max;             // Umbral superior de plausibilidad de delta
}

parameters {
  vector[P]               beta;      // Efectos fijos: intercepto, zona, luminosidad, cobertura urbana
  vector[Q]                gamma;     // Efectos aleatorios por region
  real<lower=0>             sigma;     // Escala de la log-normal (sobredispersion)
  real<lower=0>             sigmaU;    // Escala de los efectos aleatorios
  vector<lower=0>[D_obs]   densidad;  // Densidad latente (personas por VPO), segmentos aptos
}

transformed parameters {
  vector[D_obs]           lp;        // Predictor lineal, segmentos aptos
  vector<lower=0>[D_obs]  lambda;    // Intensidad Poisson, segmentos aptos

  lp = X_obs * beta + Z_obs * gamma;

  for (d in 1:D_obs)
    lambda[d] = fmin(fmax(1, densidad[d] * V_obs[d]), 10000000);
}

model {
  // Distribuciones previas (Seccion "Distribuciones previas", Capitulo 6)
  beta   ~ normal(0, 100);
  gamma  ~ normal(0, sigmaU);
  sigma  ~ cauchy(0, 2.5);
  sigmaU ~ normal(0, 10);

  // Verosimilitud (Ecuacion de verosimilitud conjunta, Capitulo 6)
  Y_obs    ~ poisson(lambda);
  densidad ~ lognormal(lp, sigma);
}

generated quantities {
  vector[D_tot]              lp_tot;       // Predictor lineal, todos los segmentos de la region
  array[D_tot] real<lower=0> densidad_hat; // Densidad predictiva, truncada en delta_max
  vector<lower=0>[D_tot]     lambda_tot;   // Intensidad Poisson predictiva
  array[D_tot] int<lower=0>  Y_tot;        // Conteo predictivo, todos los segmentos de la region
  vector[D_obs]               log_lik;      // Log-verosimilitud puntual (diagnostico LOO/WAIC)

  lp_tot = X_tot * beta + Z_tot * gamma;

  for (d in 1:D_tot) {
    densidad_hat[d] = fmin(lognormal_rng(lp_tot[d], sigma), delta_max);
    lambda_tot[d]   = fmin(fmax(1, densidad_hat[d] * V_tot[d]), 10000000);
    Y_tot[d]        = poisson_rng(lambda_tot[d]);
  }

  for (d in 1:D_obs)
    log_lik[d] = poisson_lpmf(Y_obs[d] | lambda[d]);
}

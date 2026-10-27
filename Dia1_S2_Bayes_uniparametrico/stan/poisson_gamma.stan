// Modelo Poisson-Gamma con exposicion: personas por vivienda ocupada
data {
  int<lower=0> n;              // numero de segmentos
  array[n] int<lower=0> y;     // personas censadas
  vector<lower=0>[n] v;        // viviendas particulares ocupadas (exposicion)
  real<lower=0> a;             // previa Gamma(a, b) para delta
  real<lower=0> b;
}
parameters {
  real<lower=0> delta;         // personas por vivienda
}
model {
  delta ~ gamma(a, b);
  y ~ poisson(delta * v);      // E(Y_i) = delta * V_i
}
generated quantities {
  // Prediccion para un segmento nuevo de 50 viviendas
  int y_nuevo = poisson_rng(delta * 50);
}

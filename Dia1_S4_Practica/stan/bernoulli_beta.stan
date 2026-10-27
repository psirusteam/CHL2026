// Modelo Bernoulli-Beta: proporcion theta de segmentos con cobertura completa
data {
  int<lower=0> n;                   // numero de segmentos
  array[n] int<lower=0, upper=1> y; // 1 = cobertura completa
  real<lower=0> a;                  // hiperparametros de la previa Beta(a, b)
  real<lower=0> b;
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  theta ~ beta(a, b);   // distribucion previa
  y ~ bernoulli(theta); // verosimilitud
}

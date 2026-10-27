// Modelo Multinomial-Dirichlet: distribucion de segmentos por cobertura
data {
  int<lower=1> K;              // numero de categorias
  array[K] int<lower=0> y;     // segmentos en cada categoria
  vector<lower=0>[K] alpha;    // previa Dirichlet(alpha)
}
parameters {
  simplex[K] theta;            // proporciones: positivas y suman 1
}
model {
  theta ~ dirichlet(alpha);
  y ~ multinomial(theta);
}

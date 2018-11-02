data {
  int N;
  real y[N];
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  y ~ normal(mu,sigma);
}
generated quantities {
  real y_rep[N];
  for (n in 1:N)
    y_rep[n] = normal_rng(mu,sigma);
}

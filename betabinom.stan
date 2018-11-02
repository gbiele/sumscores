data {
  int N;
  int y[N];
  int J;
}
parameters {
  real<lower=0, upper=1> p;
  real<lower=0, upper=1> theta;
}
transformed parameters {
  real<lower = 0> phi = (1-theta)/theta;
}
model {
  y ~ beta_binomial(J, p * phi, (1-p) * phi);
}
generated quantities {
  real y_rep[N];
  for (n in 1:N)
    y_rep[n] = beta_binomial_rng(J, p * phi, (1-p) * phi);
}

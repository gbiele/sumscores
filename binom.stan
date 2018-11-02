data {
  int N;
  int y[N];
  int J;
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  y ~ binomial(J,theta);
}
generated quantities {
  real y_rep[N];
  for (n in 1:N)
    y_rep[n] = binomial_rng(J,theta);
}

data {
  int K;
  int N;
  int D;
  int y[N];
  matrix[N,D] x;
}

transformed data{
  matrix[D,N] xt;

  xt = x';
}

parameters {
  matrix[K-1,D] beta_raw;

}

transformed parameters{
  matrix[K, D] beta;

  for (d in 1:D) {
    beta[1,d] = -sum(beta_raw[,d]);
    beta[2:K,d] = beta_raw[,d];
  }
}

model {
  matrix[K,N] L;

  L = beta * xt;
  
  to_vector(beta_raw) ~ normal(0, 10);

  for (n in 1:N)
    y[n] ~ categorical_logit(L[,n]);
}
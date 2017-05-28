data {
  int K;
  int N;
  int D;
  int y[N];
  matrix[N,D] x;
}

transformed data{
  vector[D] zeros;
  
  zeros = rep_vector(0, D);
}

parameters {
  matrix[D,K-1] beta_raw;
}

transformed parameters{
  matrix[D, K] beta;
  beta = append_col(zeros, beta_raw);
}

model {
  matrix[N, K] L;

  L = x * beta;
  to_vector(beta_raw) ~ normal(0, 10);
  
  for (n in 1:N)
    y[n] ~ categorical_logit(to_vector(L[n]));
}
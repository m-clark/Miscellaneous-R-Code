data {
  int<lower=1> K;                   # number of classes
  int<lower=1> N;                   # nrow of x
  int<lower=1> D;                   # ncol of x
  int<lower=1, upper=K> y[N];       # target as integer
  matrix[N,D] x;
}

transformed data{
  matrix[D,N] xt;
  row_vector[D] zeros;
  
  xt = x';
  zeros = rep_row_vector(0, D);
}

parameters {
  matrix[K-1,D] beta_raw;
}

transformed parameters{
  matrix[K, D] beta;
  beta = append_row(zeros, beta_raw);
}

model {
  matrix[K,N] L;

  L = beta * xt;
  to_vector(beta_raw) ~ normal(0, 10);
  
  for (n in 1:N)
    y[n] ~ categorical_logit(L[,n]);
}
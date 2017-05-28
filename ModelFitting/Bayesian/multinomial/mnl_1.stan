data {
  int<lower=1> K;                   # number of classes
  int<lower=1> N;                   # nrow of x
  int<lower=1> D;                   # ncol of x
  int<lower=1, upper=K> y[N];       # target as integer
  vector[D] x[N];                   # array of D
}

transformed data {
  row_vector[D] zeros;              # create reference level coefs of zero
  zeros = rep_row_vector(0, D);
}

parameters {
  matrix[K-1,D] beta_raw;           # estimated coefs
}

transformed parameters{
  matrix[K, D] beta;
  beta = append_row(zeros, beta_raw);
}

model {
  # prior
  to_vector(beta_raw) ~ normal(0, 10);
  
  # likelihood
  for (n in 1:N)
    y[n] ~ categorical_logit(beta * x[n]);
}

data {
  int K;
  int N;
  int D;
  int y[N];
  vector[D] x[N];
}

parameters {
  matrix[K,D] beta;
}

model {
  to_vector(beta) ~ normal(0, 10);
  
  for (n in 1:N)
    # y[n] ~ categorical(softmax(beta * x[n]));
  y[n] ~ categorical_logit(beta * x[n]);
}

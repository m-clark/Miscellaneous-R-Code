data {
  int N;
  int J;
  int Y[N,J];
}

transformed data{
  
}

parameters {
  vector[J] difficulty;
  real<lower=0> discrim;
  vector[N] Z;
}

transformed parameters {
  
}

model {
  matrix[N, J] lmat;
  
  # priors
  Z ~ normal(0, 1);
  discrim ~ student_t(3, 0, 5);
  difficulty ~ student_t(3, 0, 5);

  
  
  for (j in 1:J){
    lmat[,j] = discrim * (Z - difficulty[j]);
  }
  
  
  // likelihood
  for (j in 1:J)  Y[,j] ~ bernoulli_logit(lmat[,j]);
  
}

generated quantities {
  
}

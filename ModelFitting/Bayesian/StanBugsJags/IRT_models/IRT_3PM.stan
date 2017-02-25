data {
  int N;
  int J;
  int Y[N,J];
}

transformed data{
  
}

parameters {
  vector[J] difficulty;
  vector<lower=0>[J] discrim;
  vector<lower=0, upper=.25>[J] guess;
  vector[N] Z;
}

transformed parameters {
  
}

model {
  matrix[N, J] pmat;
  
  # priors
  Z ~ normal(0, 1);
  discrim ~ student_t(3, 0, 5);
  difficulty ~ student_t(3, 0, 5);
  guess ~ beta(1, 19);

  
  for (j in 1:J){
    pmat[,j] = guess[j] + (1 - guess[j]) * inv_logit(discrim[j] * (Z - difficulty[j]));
  }

  
  // likelihood
  for (j in 1:J)  Y[,j] ~ bernoulli(pmat[,j]);
  
}

generated quantities {
  
}

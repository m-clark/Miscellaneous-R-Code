data {
  int<lower=1> N;
  vector[N] X;
  vector[N] Med;
  vector[N] y;
  int<lower=1> J;                       # number of groups
  int<lower=1,upper=J> Group[N];        # Groups
}

parameters{
  real alpha_xm;                        # mediator model reg parameters
  real beta_xm;
  real<lower=0> sigma_xm;

  real alpha;                           # main model reg parameters
  real beta1;
  real beta2;
  real<lower=0> sigma_y;

  real<lower=0> sigma_alpha_xm;         # sd for ranef
  real<lower=0> sigma_alpha;
  
  vector[J] gammaAlpha_xm;
  vector[J] gammaAlpha;
}

model {
  vector[N] mu_y;                       # linear predictors for Y and M
  vector[N] mu_Med;
  
  # priors
  alpha_xm ~ normal(0, 10);   
  beta_xm ~ normal(0, 10);
  sigma_xm ~ cauchy(0,5);

  alpha ~ normal(0, 10);      
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  sigma_y ~ cauchy(0,5);

  sigma_alpha_xm  ~ cauchy(0,5);
  sigma_alpha ~ cauchy(0,5);

  for (j in 1:J){
    gammaAlpha_xm[j] ~ normal(0, sigma_alpha_xm);
    gammaAlpha[j]    ~ normal(0, sigma_alpha);
  }

  for (n in 1:N){
    mu_Med[n] <- alpha_xm + gammaAlpha_xm[Group[n]] + beta_xm*X[n];
    mu_y[n] <- alpha + gammaAlpha[Group[n]] + beta1*X[n] + beta2*Med[n] ;
  }


  # sampling
  Med ~ normal(mu_Med, sigma_xm);
  y ~ normal(mu_y, sigma_y);
}

generated quantities{
  real indEffect;

  # since no random slope for beta2
  indEffect <- beta_xm *beta2;
}
data {                                    // data setup
  int<lower=1> N;                         // sample size
  int<lower=1> I;                         // number of subjects
  vector<lower=0>[N] RT;                  // Response: reaction time
  vector<lower=0>[N] Days;                // Days in study
  int<lower=1,upper=I> Subject[N];        // Subject
}

transformed data {
  real IntBase;
  real RTsd;
  
  IntBase <- mean(RT);                    // Intercept starting point
  RTsd <- sd(RT);
}

parameters {
  real Intercept01;                       // fixed effects
  real beta01;
  vector<lower=0>[2] sigma_u;             // sd for ints and slopes
  real<lower=0> sigma_y;                  // residual sd
  vector[2] gamma[I];                     // individual effects
  cholesky_factor_corr[2] Omega_chol;     // correlation matrix for random intercepts and slopes (chol decomp)
}

transformed parameters {
  vector[I] gammaIntercept;               // individual effects (named)
  vector[I] gammaDays;
  real Intercept;
  real beta;

  Intercept <- IntBase + Intercept01 * RTsd;
  beta <- beta01 * 10;

  for (i in 1:I){  
    gammaIntercept[i]  <- gamma[i,1];
    gammaDays[i] <- gamma[i,2];
  }

} 

model {
  matrix[2,2] D;
  matrix[2,2] DC;
  vector[N] yhat;                         // Linear predictor
  vector[2] mu;                           // vector of Intercept and beta

  D <- diag_matrix(sigma_u);
  mu[1] <- Intercept;
  mu[2] <- beta;

  // priors
  Intercept01 ~ normal(0, 1);             // example of weakly informative priors;
  beta01 ~ normal(0, 1);                  // remove to essentially duplicate lme4 via improper prior

  Omega_chol ~  lkj_corr_cholesky(2.0); 

  sigma_u ~ cauchy(0, 2.5);               // prior for RE scale
  sigma_y ~ cauchy(0, 2.5);               // prior for residual scale

  DC <- D * Omega_chol;

  for (i in 1:I)                          // loop for Subject random effects
    gamma[i] ~ multi_normal_cholesky(mu, DC);

  // likelihood
  for (n in 1:N)                          
    yhat[n] <- gammaIntercept[Subject[n]] + gammaDays[Subject[n]] * Days[n];

  RT ~ normal(yhat, sigma_y);
}

generated quantities {
  matrix[2,2] Omega;                      // correlation of RE
  
  Omega <- tcrossprod(Omega_chol);
}
//see rstan_MixedModelBetaRegression.R
data {                              
  int<lower=0> N;                   // number of observations
  int<lower=1> L;                   // number of batches
  vector[N] yield;                  // response
  int<lower=1,upper=L> id[N];       // batch
  vector[N] temp;                   // temperature
}

transformed data {
  vector[N] tempCen;
  tempCen <- temp - mean(temp);     // centered explanatory variable
}

parameters {
  real Intercept;                   // "fixed" effects
  real betaTemp;

  real<lower=0> phi;                // dispersion parameter

  real<lower=0> sd_int;             // sd for ints
  real<lower=0> sd_beta;            // sd for temp

  vector[L] gammaIntercept;         // individual effects
  vector[L] gammaTemp;              // individual effects
}

transformed parameters{
  vector<lower=0>[N] A;             // parameter for beta distn
  vector<lower=0>[N] B;             // parameter for beta distn
  vector<lower=0,upper=1>[N] yhat;  // transformed linear predictor
  vector[L] IntRE;
  vector[L] SlopeRE;

  for (l in 1:L){
    IntRE[l] <- gammaIntercept[l]*sd_int;
    SlopeRE[l] <- gammaTemp[l]*sd_beta ;
  }
  
  // model calculations
  for(n in 1:N) {
    yhat[n] <- inv_logit((IntRE[id[n]] + Intercept) + (SlopeRE[id[n]] + betaTemp) * tempCen[n]);   
  }
  
  A <- yhat * phi;           
  B <- (1.0-yhat) * phi;     
}

model {
  // priors
  Intercept ~ normal(0, 10);
  betaTemp ~ normal(0, 1);          
  
  sd_int ~ cauchy(0, 2.5);        
  sd_beta ~ cauchy(0, 2.5);        
  phi ~ cauchy(0, 5);         
                                    // matt trick used for following; 
                                    // else slower and convergence issues
  gammaIntercept ~ normal(0, 1);    // random intercepts for each batch
  gammaTemp ~ normal(0, 1);         // random slopes for each batch

  // likelihood
  yield ~ beta(A, B);
}
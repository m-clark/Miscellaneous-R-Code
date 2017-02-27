#--------------------------------------------------------------------------------#
# Same as rstan_MixedModelSleepstudy.R but now estimating the correlation of the #
# random effects.  Part of this code was based on that seen on this thread       #
# https://groups.google.com/d/msg/stan-users/pdfignYQcas/BL0LPbGA2eMJ            #
#--------------------------------------------------------------------------------#

#############
### Setup ###
#############

### Data ###
### Create a model for later comparison ###
library(lme4)
data(sleepstudy)
# ?sleepstudy

mod_lme = lmer(Reaction~Days+(Days|Subject), sleepstudy)

dat = list(N=nrow(sleepstudy), I=length(unique(sleepstudy$Subject)), 
           Subject=as.numeric(sleepstudy$Subject), Days = sleepstudy$Days, 
           RT=sleepstudy$Reaction)

### Stan code ###  
stanmodelcode = '
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
  
  IntBase = mean(RT);                    // Intercept starting point
  RTsd = sd(RT);
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

  Intercept = IntBase + Intercept01 * RTsd;
  beta = beta01 * 10;

  for (i in 1:I){  
    gammaIntercept[i]  = gamma[i,1];
    gammaDays[i] = gamma[i,2];
  }

} 

model {
  matrix[2,2] D;
  matrix[2,2] DC;
  vector[N] yhat;                         // Linear predictor
  vector[2] mu;                           // vector of Intercept and beta

  D = diag_matrix(sigma_u);
  mu[1] = Intercept;
  mu[2] = beta;

  // priors
  Intercept01 ~ normal(0, 1);             // example of weakly informative priors;
  beta01 ~ normal(0, 1);                  // remove to essentially duplicate lme4 via improper prior

  Omega_chol ~  lkj_corr_cholesky(2.0); 

  sigma_u ~ cauchy(0, 2.5);               // prior for RE scale
  sigma_y ~ cauchy(0, 2.5);               // prior for residual scale

  DC = D * Omega_chol;

  for (i in 1:I)                          // loop for Subject random effects
    gamma[i] ~ multi_normal_cholesky(mu, DC);

  // likelihood
  for (n in 1:N)                          
    yhat[n] = gammaIntercept[Subject[n]] + gammaDays[Subject[n]] * Days[n];

  RT ~ normal(yhat, sigma_y);
}

generated quantities {
  matrix[2,2] Omega;                      // correlation of RE
  
  Omega = tcrossprod(Omega_chol);
}
'

#############
### Model ###
#############

library(rstan)
fit = stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 2000, warmup=200, thin=1, chains = 2,
            verbose = F)

print(fit, digits_summary=3, pars=c('Intercept', 'beta','sigma_y', 'sigma_u', 'Omega[1,2]'),
      probs = c(.025, .5, .975))

### Compare
mod_lme

print(fit, digits_summary=3, pars=c('gammaIntercept', 'gammaDays'))

### Diagnostic plots
shinystan::launch_shinystan(fit)


###############################
### A parallelized approach ###
###############################
iter = 12000
wu = 2000
thin = 10
chains = 4

fit2 = stan(model_code=stanmodelcode, model_name="mixedreg",
            fit=fit, data=dat, iter=iter, warmup=wu, 
            thin=thin, cores=4)


# some diagnostics
samplerpar = get_sampler_params(fit2)[[1]]
summary(samplerpar)


print(fit2, pars= c('Intercept', 'beta','sigma_y', 'sigma_u', 'Omega[1,2]', 'lp__'), digits=3, 
      probs = c(.025, .5, 0.975))

# Compare again
mod_lme

# Diagnostics
shinystan::launch_shinystan(fit2)


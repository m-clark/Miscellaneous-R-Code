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

parameters {
  real Intercept;                         // fixed effects
  real beta;
  vector<lower=0>[2] sigma_u;             // sd for ints and slopes
  real<lower=0> sigma_y;                  // residual sd
  vector[2] gamma[I];                     // individual effects
  corr_matrix[2] Omega;                   // correlation matrix for random intercepts and slopes
}

transformed parameters {
  vector[I] gammaIntercept;               // individual effects (named)
  vector[I] gammaDays;
  vector[N] yhat; 
  matrix[2,2] D;
  vector[2] mu;
  
  D <- diag_matrix(sigma_u);

  mu[1] <- Intercept * 100;               // see "optimizing stan code" section of manual
  mu[2] <- beta * 10;

  for (i in 1:I){  
    gammaIntercept[i]  <- gamma[i,1];
    gammaDays[i] <- gamma[i,2];
  }

  for (n in 1:N)                          // Linear predictor
    yhat[n] <- gammaIntercept[Subject[n]] + gammaDays[Subject[n]] * Days[n];
} 

model {
  matrix[2,2] C;                          // for cholesky decomposition of corr matrix
  matrix[2,2] DC;

  // priors
  Intercept ~ normal(0, 1);               // example of weakly informative priors;
  beta ~ normal(0, 1);                    // remove to essentially duplicate lme4 via improper prior

  Omega ~  lkj_corr(1.0); 

  sigma_u ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);

  C <- cholesky_decompose(Omega);
  DC <- D * C;

  for (i in 1:I)                          // loop for Subject random effects
    gamma[i] ~ multi_normal_cholesky(mu, DC);

  // likelihood
  RT ~ normal(yhat, sigma_y);
}
'

#############
### Model ###
#############

library(rstan)
fit = stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 2000, warmup=200, thin=1, chains = 2,
            verbose = F)

print(fit, digits_summary=3, pars=c('mu','sigma_y', 'sigma_u', 'Omega[1,2]'),
      probs = c(0, .025, .5, .975, 1))

### Compare
mod_lme

print(fit, digits_summary=3, pars=c('gammaIntercept', 'gammaDays'))

### Diagnostic plots
traceplot(fit, pars=c('mu','sigma_y', 'sigma_u', 'Omega[1,2]'))


###############################
### A parallelized approach ###
###############################
library(parallel)
cl = makeCluster(3)
clusterEvalQ(cl, library(rstan))

clusterExport(cl, c('stanmodelcode', 'dat', 'fit')) 

p = proc.time()
parfit = parSapply(cl, 1:3, function(i) stan(model_code = stanmodelcode, model_name = "mixedreg", #init=0,
                                              fit = fit, # if using the same model code, this will use the previous compilation
                                              data = dat, iter = 120000, warmup=20000, thin=10, chains = 1, chain_id=i,
                                              control=list(stepsize=.01, stepsize_jitter=.5, max_treedepth=20),
                                              verbose = T), 
                    simplify=F) 

proc.time() - p

stopCluster(cl)

# combine the chains
fit2 = sflist2stanfit(parfit)

# examine some diagnostics
ainfo = get_adaptation_info(fit2)
cat(ainfo[[1]])
samplerpar = get_sampler_params(fit2)[[1]]
summary(samplerpar)


print(fit2, pars= c('mu','sigma_y', 'sigma_u', 'Omega[1,2]', 'lp__'), digits=3, 
      probs = c(.01, .025, .05, .5, .95, 0.975, .99))

# Compare again
mod_lme

# Diagnostics
traceplot(fit2, inc_warmup=F, pars=c('mu','sigma_y', 'sigma_u', 'Omega[1,2]', 'lp__'))

pairs(fit2, pars=c('sigma_u'))
pairs(fit2, pars=c('mu'))
# plot(fit2)

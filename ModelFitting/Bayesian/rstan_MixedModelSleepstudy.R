#---------------------------------------------------------------------------------#
# A mixed model via stan/rstan with comparison to lme4 output.  In this model     #
# there is a random intercept as well as a random slope for the single predictor. #
#---------------------------------------------------------------------------------#

#############
### Setup ###
#############

### Data ###
library(lme4)
data(sleepstudy)
# ?sleepstudy

# Create a model for later comparison
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
  real<lower=0> sd_int;                   // sd for ints and slopes
  real<lower=0> sd_beta;
  real<lower=0> sigma_y;                  // residual sd

  vector[I] gammaIntercept;               // individual effects
  vector[I] gammaDays;
}

transformed parameters {
  vector[N] yhat; 

  for (n in 1:N)                         // Linear predictor
    yhat[n] = gammaIntercept[Subject[n]] + gammaDays[Subject[n]] * Days[n];
} 

model {
  // priors
  Intercept ~ normal(0, 100);            // example of weakly informative priors (and ignoring Matt trick for now);
  beta ~ normal(0, 100);                 // remove to essentially duplicate lme4 via improper prior

  gammaIntercept ~ normal(Intercept, sd_int);
  gammaDays ~ normal(beta, sd_beta);
  
  sd_int ~ cauchy(0, 2.5);
  sd_beta ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);

  // likelihood
  RT ~ normal(yhat, sigma_y);
}
'
#############
### Model ###
#############

library(rstan)
fit = stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 22000, warmup=2000, thin=20, chains = 4,
            verbose = F) 

### Summarize
print(fit, digits_summary=3, pars=c('Intercept','beta','sigma_y', 'sd_int', 'sd_beta'),
      probs = c(0, .025, .5, .975, 1))

### Compare
mod_lme

print(fit, digits_summary=3, pars=c('gammaIntercept', 'gammaDays'))

### Diagnostic plots
shinystan::launch_shinystan(fit)

###############################
### A parallelized approach ###
###############################

fit2 = stan(model_code = stanmodelcode, model_name = "mixedreg", #init=0,
            fit = fit, # if using the same model code, this will use the previous compilation
            data = dat, iter = 120000, warmup=20000, thin=10, cores=3,
            verbose = T)


# examine some diagnostics
ainfo = get_adaptation_info(fit2)
cat(ainfo[[1]])
samplerpar = get_sampler_params(fit2)[[1]]
summary(samplerpar)


print(fit2, pars= c('Intercept','beta','sigma_y', 'sd_int', 'sd_beta','lp__'), digits=3, 
      probs = c(.01, .025, .05, .5, .95, 0.975, .99))

# Compare again
mod_lme

# diagnostics
shinystan::launch_shinystan(fit2)

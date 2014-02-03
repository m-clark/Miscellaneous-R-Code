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

# Make a meaningful intercept
Days2 = sleepstudy$Days-1

# Create a model for later comparison
mod_lme = lmer(Reaction~Days2+(Days2|Subject), sleepstudy)

X = model.matrix(mod_lme)
y = sleepstudy[,1]
Subject = sleepstudy$Subject

dat = list(N=nrow(sleepstudy), I=length(unique(sleepstudy$Subject)), 
           Subject=as.numeric(sleepstudy$Subject), Days = Days2, 
           RT=sleepstudy$Reaction)




### Stan code ###
stanmodelcode <-'
data {                                    // data setup
  int<lower=1> N;                         // sample size
  int<lower=1> I;                         // number of subjects
  vector[N] RT;                           // Response- reaction time
  vector[N] Days;                         // Days in study
  int<lower=1,upper=I> Subject[N];        // Subject
}

parameters {
  real Int;                               // fixed effects
  real beta;
  real<lower=0> sd_int;                   // sd for ints and slopes
  real<lower=0> sd_beta;
  real<lower=0> sigma_y;                  // residual sd

  vector[I] gammaInt;                     // individual effects
  vector[I] gammaDays;
}

transformed parameters {
  vector[N] yhat; 

  for (n in 1:N)                         // Linear predictor
    yhat[n] <- gammaInt[Subject[n]] + gammaDays[Subject[n]] * Days[n];
} 

model {
  // priors
  Int ~ normal(0, 100);                  // example of weakly informative priors (and ignoring Matt trick for now);
  beta ~ normal(0, 100);                 // remove to essentially duplicate lme4 via improper prior

  gammaInt ~ normal(Int, sd_int);
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
fit <- stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 22000, warmup=2000, thin=20, chains = 4,
            verbose = F) 

### Summarize
print(fit, digits_summary=3, pars=c('Int','beta','sigma_y', 'sd_int', 'sd_beta'),
      probs = c(0, .025, .5, .975, 1))

### Compare
mod_lme

print(fit, digits_summary=3, pars=c('gammaInt', 'gammaDays'))

### Diagnostic plots
traceplot(fit, pars=c('Int','beta','sigma_y', 'sd_int', 'sd_beta'))
traceplot(fit, pars=c('Int','beta','sigma_y', 'sd_int', 'sd_beta'), inc_warmup=F)

###############################
### A parallelized approach ###
###############################
library(parallel)
cl = makeCluster(3)
clusterEvalQ(cl, library(rstan))

clusterExport(cl, c('stanmodelcode', 'dat', 'fit')) 

p = proc.time()
parfit <- parSapply(cl, 1:3, function(i) stan(model_code = stanmodelcode, model_name = "mixedreg", #init=0,
                                                   fit = fit, # if using the same model code, this will use the previous compilation
                                                   data = dat, iter = 120000, warmup=20000, thin=10, chains = 1, chain_id=i,
                                                   verbose = T), 
                         simplify=F) 

proc.time() - p

stopCluster(cl)

# combine the chains
fit2 = sflist2stanfit(parfit)

# examine some diagnostics
ainfo <- get_adaptation_info(fit2)
cat(ainfo[[1]])
samplerpar = get_sampler_params(fit2)[[1]]
summary(samplerpar)


print(fit2, pars= c('Int','beta','sigma_y', 'sd_int', 'sd_beta','lp__'), digits=3, 
      probs = c(.01, .025, .05, .5, .95, 0.975, .99))

# Compare again
mod_lme

library(scales)
traceplot(fit2, inc_warmup=F, pars=c('Int','beta','sigma_y', 'sd_int', 'sd_beta', 'lp__'))
pairs(fit2, pars=c('Int','beta'))

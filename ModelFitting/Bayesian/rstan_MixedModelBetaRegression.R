####################################################################################
### The following is an attempt to use Stan/rstan for a mixed model assuming a   ###
### beta distribution for the response. See the helpfile for GasolineYield for   ###
### information on the data.  A similar example in JAGS can be found in          ###
### jags_MixedModelBetaRegression.R. In general you'd probably want tweak the    ###
### settings for both some more for optimal results, but for demonstration the   ###
### following will produce adequate results.                                     ###
####################################################################################

####################
### Get the data ###
####################
data(GasolineYield, package='betareg')
hist(GasolineYield$yield)


# for a quick comparison we might use the mgcv package; note that mgcv seems to 
# be the closest one can get to an 'out-of-the-package' mixed model for a beta 
# distributed response in R (maybe try gamlss, but its results seem problematic).
# However, for random slopes, e.g. with gamm4/gamm, you no longer have access to
# the beta distribution due to underlying use of lme4/nlme. The addition of
# random slopes doesn't seem to have too much effect here, so it can still be a
# useful comparison to see if it and Stan are in the same ball park as far as
# the random intercepts go; otherwise take out the random slopes of the Stan
# model for a direct comparison, and you'll get practically identical estimates.
# As mentioned above, a JAGS version of the model is available for direct
# comparison.

library(mgcv)
modgam = gam(yield ~ scale(temp, scale=F) + s(batch, bs='re'), data=GasolineYield, family=betar(link="logit"))
# summary(modgam)  

#######################
### Stan model code ###
#######################
stanmodelcode <-'
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
'


################
### Test run ###
################

### Data set up
library(rstan)
modelMatrix = data.frame(1, GasolineYield[,'temp', drop=F])

dat = list(N = nrow(GasolineYield), yield=GasolineYield$yield, 
           temp=GasolineYield$temp, id=as.numeric(GasolineYield$batch), L=10) # 10 = number of batches


### Run the model and check diagnostics
test <- stan(model_code = stanmodelcode, model_name = "example", 
             init=0,  # initializing unbounded at zero or setting the range for random init values proved useful
             # init_r=.01,  
             data = dat, iter = 2000, warmup=200, thin=1, chains = 2, 
             pars = c('Intercept','betaTemp','IntRE', 'SlopeRE', 'phi', 'sd_int', 'sd_beta'), verbose = F, refresh=2000/4) 

# diagnostics
traceplot(test, inc=T)
ainfo <- get_adaptation_info(test)
cat(ainfo[[1]])
samplerpar = get_sampler_params(test)[[1]]
summary(samplerpar)
pairs(test, pars=c('Intercept', 'betaTemp'))

# summary
print(test, digits_summary=4)



################
### Full run ###
################

### set iterations etc.
iters = 22000
wu = 2000
thin = 20
nchains = 4


### setup and run in parallel
library(parallel)
cl = makeCluster(4)

clusterEvalQ(cl, library(rstan))
clusterExport(cl, c('stanmodelcode', 'dat', 'test', 'iters', 'wu', 'thin')) 

p = proc.time()
parfit = parSapply(cl, 1:nchains, function(i) stan(model_code = stanmodelcode, model_name = "mixedreg", 
                                                   init=0,       # see comment in test run
                                                   # init_r=.01, 
                                                   fit = test, 
                                                   data = dat, iter = iters, warmup=wu, thin=thin, chains = 1, chain_id=i,
                                                   verbose = T, refresh=iters/4), 
                   simplify=F) 

proc.time() - p

stopCluster(cl)


### combine the chains
fit2 = sflist2stanfit(parfit)


### examine some diagnostics
ainfo = get_adaptation_info(fit2)
cat(ainfo[[1]])
samplerpar = get_sampler_params(fit2)[[1]]
summary(samplerpar)


### examine model
print(fit2, pars = c('Intercept','betaTemp','IntRE', 'SlopeRE', 'phi', 'sd_int', 'sd_beta', 'lp__'), 
      digits=4, probs = c(.025, .5, 0.975))

traceplot(fit2, inc_warmup=F, pars = c('Intercept','betaTemp','IntRE', 'SlopeRE', 'phi', 'sd_int', 'sd_beta', 'lp__'))


### Compare to mgcv
summary(modgam)
gam.vcomp(modgam)
gamIntRE = modgam$coeff[grep('batch', names(modgam$coef))]
stanIntRE = colMeans(extract(fit2, par='IntRE')[[1]])
cbind(gamIntRE, stanIntRE)



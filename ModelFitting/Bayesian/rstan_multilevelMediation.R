# The following demonstrates an indirect effect in a multilevel situation.  It 
# is based on Yuan & MacKinnon 2009, which provides some Bugs code.  In what
# follows we essentially have two models, one where the 'mediator' is the
# response; the other regards the primary response of interest (noted y).  They will be
# referred to with Med or Main respectively.  I also don't follow the notation
# in the article as I don't find it as clear.
set.seed(8675309)

##################
### Data Setup ###
##################
### Parameters
## the two main models are expressed conceptually as
# Mediator ~ alphaMed + betaMed*X
# y ~ alphaMain + beta1Main*X + beta2Main*Mediator

# additionally there will be random effects for a grouping variable for each
# coefficient, i.e. random intercepts and slopes

library(MASS) # for mvrnorm

## random effects for mediator model
# create cov matrix of RE etc. with no covariance between model random effects
# covmat_RE = matrix(c(1,-.15,0,0,0,
#                       -.15,.4,0,0,0,
#                       0,0,1,-.1,.15,
#                       0,0,-.1,.3,0,
#                       0,0,.15,0,.2), nrow=5, byrow = T)

# or with slight cov added to indirect coefficient RE; both matrices are pos def
covmat_RE = matrix(c(1,-.15,0,0,0,
                     -.15,.64,0,0,-.1,
                     0,0,1,-.1,.15,
                     0,0,-.1,.49,0,
                     0,-.1,.15,0,.25), nrow=5, byrow = T)

# inspect
covmat_RE

# inspect as correlation
cov2cor(covmat_RE)

# simulate
re = mvrnorm(50, mu=rep(0,5), Sigma=covmat_RE, empirical = T)

# random effects for mediator model
ranef_alphaMed = rep(re[,1], e=10)
ranef_betaMed = rep(re[,2], e=10)

# random effects for main model                                                 
ranef_alphaMain = rep(re[,3], e=10)
ranef_beta1Main = rep(re[,4], e=10)
ranef_beta2Main = rep(re[,5], e=10)

## fixed effects
alphaMed = 2
betaMed = .2

alphaMain = 1
beta1Main = .3
beta2Main = -.2

# residual variance
resid_Med = mvrnorm(500, 0, .75^2, empirical=T)
resid_Main = mvrnorm(500, 0, .5^2, empirical=T)


# Collect parameters for later comparison
params = c(alphaMed=alphaMed, betaMed=betaMed, sigmaMed=sd(resid_Med), 
           alphaMain=alphaMain, beta1Main=beta1Main, beta2Main=beta2Main, sigma_y=sd(resid_Main), 
           alphaMed_sd=sqrt(diag(covmat_RE)[1]), betaMed_sd=sqrt(diag(covmat_RE)[2]), 
           alpha_sd=sqrt(diag(covmat_RE)[3]), beta1_sd=sqrt(diag(covmat_RE)[4]), beta2_sd=sqrt(diag(covmat_RE)[5])
           )

ranefs =  cbind(gammaAlphaMed=unique(ranef_alphaMed), gammaBetaMed=unique(ranef_betaMed), 
                gammaAlpha=unique(ranef_alphaMain), gammaBeta1=unique(ranef_beta1Main), gammaBeta2=unique(ranef_beta2Main))


### Create Data
X = rnorm(500, sd=2)
Med = (alphaMed + ranef_alphaMed) + (betaMed+ranef_betaMed)*X + resid_Med[,1]
y = (alphaMain + ranef_alphaMain) + (beta1Main+ranef_beta1Main)*X + (beta2Main + ranef_beta2Main)*Med + resid_Main[,1]
group = rep(1:50, e=10)


### A piecemeal lme for comparison; can't directly estimate mediated effect, and it
### won't pick up on correlation of random effects between models
library(lme4)
modMed = lmer(Med ~ X + (1+X|group))
summary(modMed)

modMain = lmer(y ~ X + Med + (1+X+Med|group))
summary(modMain)

# should equal the naive estimate in the following code
medLme = fixef(modMed)[2]*fixef(modMain)[3]



#################
### Stan Code ###
#################

# In the following, the cholesky decomposition of the RE covariance matrix is
# used for efficiency.  As a rough guide, the default data where N=500 took
# about 5 min to run for the main model with iter=12000 and warmup of 2000.

modelStan <- "
data {
  int<lower=1> N;                                # Sample size
  vector[N] X;                                   # Explanatory variable
  vector[N] Med;                                 # Mediator
  vector[N] y;                                   # response
  int<lower=1> J;                                # number of groups
  int<lower=1,upper=J> Group[N];                 # Groups
}

parameters{
  real alphaMed;                                 # mediator model reg parameters and related
  real betaMed;
  real<lower=0> sigma_alphaMed;
  real<lower=0> sigma_betaMed;
  real<lower=0> sigmaMed;

  real alphaMain;                                # main model reg parameters and related
  real beta1Main;
  real beta2Main;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_beta2;
  real<lower=0> sigma_y;

  cholesky_factor_corr[5] Omega_chol;            # chol decomp of corr matrix for random effects

  vector<lower=0>[5] sigma_ranef;                # sd for random effects

  matrix[J,5] gamma;                             # random effects
}

transformed parameters{
  vector[J] gammaAlphaMed;
  vector[J] gammaBetaMed;

  vector[J] gammaAlpha;
  vector[J] gammaBeta1;
  vector[J] gammaBeta2;

  for (j in 1:J){
    gammaAlphaMed[j] <- gamma[j,1];
    gammaBetaMed[j] <- gamma[j,2];
    gammaAlpha[j] <- gamma[j,3];
    gammaBeta1[j] <- gamma[j,4];
    gammaBeta2[j] <- gamma[j,5];
  }
}

model {
  vector[N] mu_y;                                # linear predictors for response and mediator
  vector[N] mu_Med;
  matrix[5,5] D;
  matrix[5,5] DC;

  ### priors
  ## mediator model
  # fixef
  sigma_alphaMed ~ cauchy(0, 1);                 # for scale params the cauchy is a little more informative here due to the nature of the data
  sigma_betaMed ~ cauchy(0, 1);
  alphaMed ~ normal(0, sigma_alphaMed);   
  betaMed ~ normal(0, sigma_betaMed);

  # residual scale
  sigmaMed ~ cauchy(0, 1);

  ## main model
  # fixef
  sigma_alpha ~ cauchy(0, 1);
  sigma_beta1 ~ cauchy(0, 1);
  sigma_beta2 ~ cauchy(0, 1);
  alphaMain ~ normal(0, sigma_alpha);      
  beta1Main ~ normal(0, sigma_beta1);
  beta2Main ~ normal(0, sigma_beta2);

  # residual scale
  sigma_y ~ cauchy(0, 1);

  ## ranef sampling via cholesky decomposition
  sigma_ranef  ~ cauchy(0, 1);
  Omega_chol ~  lkj_corr_cholesky(2.0);

  D  <- diag_matrix(sigma_ranef);
  DC <- D * Omega_chol;
  
  for (j in 1:J)                                 # loop for Group random effects
    gamma[j] ~ multi_normal_cholesky(rep_vector(0, 5), DC);

  ## Linear predictors
  for (n in 1:N){
    mu_Med[n] <- alphaMed + gammaAlphaMed[Group[n]] + (betaMed + gammaBetaMed[Group[n]])*X[n];
    mu_y[n]   <- alphaMain + gammaAlpha[Group[n]] + (beta1Main+gammaBeta1[Group[n]])*X[n] + (beta2Main+gammaBeta2[Group[n]])*Med[n] ;
  }
  
  
  ### sampling for primary models
  Med ~ normal(mu_Med, sigmaMed);
  y ~ normal(mu_y, sigma_y);
}

generated quantities{
  real naiveIndEffect;
  real avgIndEffect;
  real totalEffect;
  matrix[5,5] covMatRE;
  
  covMatRE <- diag_matrix(sigma_ranef) * tcrossprod(Omega_chol) * diag_matrix(sigma_ranef);

  naiveIndEffect <- betaMed*beta2Main;
  avgIndEffect <- betaMed*beta2Main + covMatRE[2,5];
  totalEffect <- avgIndEffect + beta1Main;
}
"


###########################
### Debugging and Setup ###
###########################

standat = list(X=X, Med=Med, y=y, Group=group, J=length(unique(group)), N=length(y))

# bug detector
library(rstan)
p = proc.time()
fit0 = stan(model_code = modelStan, data = standat, iter = 400, warmup=200, 
            thin=1, chains = 1, verbose = F)
proc.time() - p
 
# print(fit0)

# traceplot(fit0)



##################
### Main Model ###
##################
### Setup
iter = 12000
wu = 2000
thin = 20
chains = 4

library(parallel)
cl = makeCluster(chains)
clusterEvalQ(cl, library(rstan))

clusterExport(cl, c('modelStan', 'standat', 'fit0', 'iter', 'wu', 'thin', 'chains')) 

### Run
p = proc.time()
parfit=parSapply(cl, 1:chains, function(i) stan(model_code=modelStan, fit=fit0, 
                                                data=standat, iter=iter, warmup=wu, 
                                                thin=thin, chains=1, chain_id=i), 
                 simplify=F) 

(proc.time() - p)[3]/60

stopCluster(cl)

# combine the chains
fit = sflist2stanfit(parfit)



#########################
### Model Exploration ### 
#########################
### Summarize model

# main parameters include fixed and random effect sd, plus those related to
# indirect effect
mainpars = c('alphaMed', 'betaMed', 'sigmaMed',
             'alphaMain', 'beta1Main', 'beta2Main', 'sigma_y',
             'sigma_ranef',
             'naiveIndEffect', 'avgIndEffect','totalEffect')

print(fit, digits=3, probs = c(.025, .5, 0.975), pars=mainpars)


### Compare parameters of interest

## Extract parameters for comparison
pars1 = get_posterior_mean(fit, pars=mainpars)[,5]
parsREcov = get_posterior_mean(fit, pars='Omega_chol')[,5] # or take 'covMatRE' from monte carlo sim
parsRE = get_posterior_mean(fit, pars=c('sigma_ranef'))[,5]


## fixed effects and re variances
cbind(params, pars1[1:12], c(modMed@beta, summary(modMed)$sigma, 
                       modMain@beta, summary(modMain)$sigma,
                       summary(modMed)$varcor$group[1,1]^.5,summary(modMed)$varcor$group[2,2]^.5,
                       summary(modMain)$varcor$group[1,1]^.5,summary(modMain)$varcor$group[2,2]^.5,
                       summary(modMain)$varcor$group[3,3]^.5)
      )

## compare covariance of random effects
covmat_RE_est = diag(parsRE) %*% tcrossprod(matrix(parsREcov, ncol = 5, byrow = T)) %*% diag(parsRE)
list(covmat_RE, round(covmat_RE_est, 2))

vcovMed = covmat_RE_est[1:2,1:2]
list(covmat_RE[1:2,1:2], round(vcovMed, 2), round(summary(modMed)$varcor$group[1:2,1:2], 2))

vcovMain = covmat_RE_est[3:5,3:5]
list(covmat_RE[3:5,3:5], round(vcovMain, 2), round(summary(modMain)$varcor$group[1:3,1:3], 2))


## compare indirect effects
c(true = betaMed*beta2Main + covmat_RE[2,5],
  est = get_posterior_mean(fit, 'avgIndEffect')[,5],
  naiveBayes = get_posterior_mean(fit, 'naiveIndEffect')[,5],
  naiveLme = medLme)



########################
### Diagnostics etc. ###
########################
samplerpar = get_sampler_params(fit)[[1]]
summary(samplerpar)

traceplot(fit, pars=mainpars, inc_warmup=F)
plot(fit, pars=mainpars)




# 3rd variant
covmat_RE = matrix(c(4,-1.5,0,0,0,
                     -1.5,2,0,0,-.5,
                     0,0,4,-1,1,
                     0,0,-1,2,0,
                     0,-.5,1,0,1), nrow=5, byrow = T)
# inspect
covmat_RE

# inspect as correlation
cov2cor(covmat_RE)

# simulate
re = mvrnorm(50, mu=rep(0,5), Sigma=covmat_RE, empirical = T)

# random effects for mediator model
ranef_alphaMed = rep(re[,1], e=10)
ranef_betaMed = rep(re[,2], e=10)

# random effects for main model                                                 
ranef_alphaMain = rep(re[,3], e=10)
ranef_beta1Main = rep(re[,4], e=10)
ranef_beta2Main = rep(re[,5], e=10)

## fixed effects
alphaMed = 4
betaMed = 2

alphaMain = 3
beta1Main = 1
beta2Main = -2

# residual variance
resid_Med = mvrnorm(500, 0, 2^2, empirical=T)
resid_Main = mvrnorm(500, 0, 1^2, empirical=T)


# Collect parameters for later comparison
params = c(alphaMed=alphaMed, betaMed=betaMed, sigmaMed=sd(resid_Med), 
           alphaMain=alphaMain, beta1Main=beta1Main, beta2Main=beta2Main, sigma_y=sd(resid_Main), 
           alphaMed_sd=sqrt(diag(covmat_RE)[1]), betaMed_sd=sqrt(diag(covmat_RE)[2]), 
           alpha_sd=sqrt(diag(covmat_RE)[3]), beta1_sd=sqrt(diag(covmat_RE)[4]), beta2_sd=sqrt(diag(covmat_RE)[5])
)

ranefs =  cbind(gammaAlphaMed=unique(ranef_alphaMed), gammaBetaMed=unique(ranef_betaMed), 
                gammaAlpha=unique(ranef_alphaMain), gammaBeta1=unique(ranef_beta1Main), gammaBeta2=unique(ranef_beta2Main))


### Create Data
X = rnorm(500, sd=2)
Med = (alphaMed + ranef_alphaMed) + (betaMed+ranef_betaMed)*X + resid_Med[,1]
y = (alphaMain + ranef_alphaMain) + (beta1Main+ranef_beta1Main)*X + (beta2Main + ranef_beta2Main)*Med + resid_Main[,1]
group = rep(1:50, e=10)


### A piecemeal lme for comparison; can't directly estimate mediated effect, and it
### won't pick up on correlation of random effects between models
library(lme4)
modMed = lmer(Med ~ X + (1+X|group))
summary(modMed)

modMain = lmer(y ~ X + Med + (1+X+Med|group))
summary(modMain)

# should equal the naive estimate in the following code
medLme = fixef(modMed)[2]*fixef(modMain)[3]

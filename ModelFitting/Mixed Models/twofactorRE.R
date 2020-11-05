###################################################################################
### An approach for two factor random effects model via maximum likelihood in R ###
### Matlab and Julia; It's based on Statistical Modeling and Computation (2014) ###
### Chapter 10, example 10.10; The data regards the breeding value of a set of  ### 
### five sires in raising pigs. Each sire is mated to a random group of dams,   ###
### with the response being the average daily weight gain in pounds of two      ###
### piglets in each litter.  See one_factor_RE.R for a one factor model, and    ###
### two_factor_RE.m and two_factor_RE.jl for the matlab and julia versions of   ###
### this example. Note that the text has a typo on the sigma2 variance estimate ###
### (value should be .0023 not .023).                                           ###
###################################################################################



#####################
### Main function ###
#####################
# the function takes the log variances eta* as input to keep positive
sfran_loglike = function(mu, eta_alpha, eta_gamma, eta) {
  # Args 
  # mu: intercept 
  # eta_alpha: random effect one 
  # eta_gamma: random effect two
  # eta: residual variance of y
  
  sigma2_alpha = exp(eta_alpha)
  sigma2_gamma = exp(eta_gamma)
  sigma2 = exp(eta)
  n = length(y)
  
  # covariance matrix of observations
  Sigma = sigma2 * diag(n) + sigma2_alpha * tcrossprod(Xalpha) + 
    sigma2_gamma * tcrossprod(Xgamma)
  
  
  # log likelihood
  ll = -n / 2 * log(2 * pi) - sum(log(diag(chol(Sigma)))) -
    .5 * t(y - mu) %*% chol2inv(chol(Sigma)) %*% (y - mu)
  
  return(-ll)
}


###################
### Data set up ###
###################
y = c(1.39,1.29,1.12,1.16,1.52,1.62,1.88,1.87,1.24,1.18,
      .95,.96,.82,.92,1.18,1.20,1.47,1.41,1.57,1.65)

# for use in lme4, but also a more conceptual respresentation of the data
d = expand.grid(sire = rep(1:5, 2), dam = 1:2)
d = data.frame(d[order(d$sire), ], y)


################################
### Starting values and test ###
################################
starts = list(
  mu = mean(y),
  eta_alpha = var(tapply(y, d$sire, mean)),
  eta_gamma = var(y) / 3,
  eta = var(y) / 3
)

Xalpha = diag(5) %x% rep(1,4)

Xgamma = diag(10) %x% rep(1,2)

# test
sfran_loglike(starts[[1]], starts[[2]], starts[[3]], starts[[4]])


#######################
### Run and compare ###
#######################
library(bbmle)

mlout = mle2(sfran_loglike, start=starts,  method='BFGS')  

### lme4 comparison
library(lme4)

lme = lmer(y ~ (1 | sire) + (1 | dam:sire), d, REML = F)

summary(mlout)
exp(coef(mlout)[-1])

summary(lme)




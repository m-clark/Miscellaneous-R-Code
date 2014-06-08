###################################################################################
### ML approach for one factor random effects model via maximum likelihood in R ###
### Matlab and Julia; It's based on Statistical Modeling and Computation (2014) ###
### Chapter 10, example 10.10; Unfortunately I did this before knowing they had ###
### both matlab and R code on their website, though the R code here is a little ###
### cleaner and has comments. The data regards crop yield from 10 randomly      ###
### selected locations and three collections at each.  See onefactorRE.m and    ###
### onefactorRE.jl for the related Matlab and Julia files, and the respective   ###
### twofactorRE.* for the associated two factor random effects examples.        ###
###################################################################################



#####################
### Main function ###
#####################
sfran_loglike = function(mu, sigma2_mu, sigma2){
  # Args are mu: intercept; sigma2_mu: variance of intercept; sigma2: residual
  # variance of y
  # I follow their (unfortunate) notation. 
  d = nrow(y)
  ni = ncol(y)
  
  # covariance matrix of observations
  Sigmai = sigma2*diag(ni) + sigma2_mu*matrix(1,ni,ni)
  
  # log likelihood
  l = rep(NA, 10)
  # iterate over the rows
  for(i in 1:d){
    l[i] = .5 * t(y[i,]-mu) %*% chol2inv(chol(Sigmai)) %*% (y[i,]-mu)  
  }
  
  ll =  -(ni*d)/2*log(2*pi) - d/2*log(det(Sigmai)) - sum(l)
  return(-ll)
}


###################
### Data set up ###
###################
y = matrix(c(22.6,20.5,20.8,
             22.6,21.2,20.5,
             17.3,16.2,16.6,
             21.4,23.7,23.2,
             20.9,22.2,22.6,
             14.5,10.5,12.3,
             20.8,19.1,21.3,
             17.4,18.6,18.6,
             25.1,24.8,24.9,
             14.9,16.3,16.6), 10, 3, byrow=T)


################################
### Starting values and test ###
################################
starts = list(mu=mean(y), sigma2_mu=var(rowMeans(y)), sigma2=mean(apply(y, 1, var)))

### test
sfran_loglike(starts[[1]], starts[[2]], starts[[3]])


#######################
### Run and compare ###
#######################

### bbmle has mle2 function for maximum likelihood estimation based on
### underlying R functions like optim. LBFGS-B is used to place lower bounds on
### the variance estimates
library(bbmle)
mlout = mle2(sfran_loglike, start=starts,  method='L-BFGS-B', 
             lower=c(mu=-Inf, sigma2_mu=0, sigma2=0), trace=T)  

### Compare
library(lme4); library(reshape2)
lme = lmer(value~1|Var1, data=melt(y), REML=F)

summary(mlout)
summary(lme); -2*logLik(lme)

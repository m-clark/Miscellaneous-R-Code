library(mvtnorm)
library(psych)


###########################
### Create the data set ###
###########################
set.seed(123)

lambda = (matrix(c(1,.5,.3,.6,0,0,0,0,
                   0,0,0,0,1,.7,.4,.5), 
                 nrow=2, byrow=T))

phi = matrix(c(1,.25,.25,1),nrow=2,byrow=T)  #correlation of factors

factors = rmvnorm(1000, mean=rep(0,2), sigma=phi, "chol")
e = rmvnorm(1000, sigma=diag(8))

y = 0 + factors%*%lambda + e

#dim(y)
describe(y)
round(cor(y),3)

#see the factor structure
library(corrplot)
corrplot(cor(y))

#fa(y,nfactors=2, rotate="oblimin") #example exploratory

#########################
### Primary Functions ###
#########################

#measurement model, covariance approach
cfa.cov <- function (parms, data) {
  S = cov(data)*((nrow(data)-1)/nrow(data))  #ML covariance div by N rather than N-1, the multiplier adjusts
  lambda = matrix(c(1,parms[1:3],rep(0,4),rep(0,4), 1,parms[4:6]), ncol=2) #loading estimates
  dist.init = parms[7:14]    #disturbances
  disturbs = diag(dist.init)
  phi.init = matrix(c(parms[15],parms[17],parms[17],parms[16]), ncol=2)  #factor cov/correlation matrix
  sigtheta = lambda%*%phi.init%*%t(lambda)+disturbs
  pq = dim(data)[2]  #in Bollen p + q (but for the purposes of thishere just p) = tr(data)
  #out = -(log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - pq )  #a reduced version; Bollen 1989 p.107
  out = ((-(nrow(data))*pq/2)*log(2*pi)) - ((nrow(data))/2)* (log(det(sigtheta)) + tr(S%*%solve(sigtheta)))  #should be same as Mplus H0 loglike
  out
}

#correlation approach for standardized output
cfa.cor = function (parms, data) {
  S = cor(data)
  lambda = matrix(c(parms[1:4],rep(0,4),rep(0,4), parms[5:8]), ncol=2)
  dist.init = parms[9:16]
  disturbs = diag(dist.init)
  phi.init = matrix(c(1,parms[17],parms[17],1), ncol=2)
  sigtheta = lambda%*%phi.init%*%t(lambda)+disturbs
  pq = dim(data)[2]
  #out = (log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - pq )
  out = ((-nrow(data)*pq/2)*log(2*pi))-(nrow(data)/2)* (log(det(sigtheta)) + tr(S%*%solve(sigtheta)))
  out
}

####################
### optimization ###
####################

par.init.cov=c(rep(0,6), rep(.9,8), rep(.5,3)) #cov
out.cov = optim(par = par.init.cov, fn=cfa.cov, data=y, method="L-BFGS-B", lower = 0, control = c(list(fnscale = -1))) #, control = c(list(fnscale = 1), parscale = rep(0.01, length(par.init)))
loads.cov= matrix(c(1,out.cov$par[1:3],rep(0,4), rep(0,4),1,out.cov$par[4:6]), ncol=2)
disturbs.cov = out.cov$par[7:14]
#round(data.frame(loads.cov,disturbs.cov),3); round(matrix(c(out.cov$par[c(15,17,17,16)]), ncol=2),3)

par.init.cor=c(rep(0,8), rep(.9,8), 0) #for cor
out.cor = optim(par = par.init.cor, fn=cfa.cor, data=y, method="L-BFGS-B", lower = 0, upper = 1, control = c(list(fnscale = -1))) #, control = c(list(fnscale = 1), parscale = rep(0.01, length(par.init)))
loads.cor= matrix(c(out.cor$par[1:4],rep(0,4), rep(0,4),out.cor$par[5:8]), ncol=2)
disturbs.cor = out.cor$par[9:16]
#round(data.frame(loads.cor,disturbs.cor),3); round(matrix(c(1,out.cor$par[c(17,17)],1), ncol=2),3)


#################################
### Gather output for summary ###
#################################

output=list(out.raw = round(data.frame(loads.cov,disturbs.cov),3), cov.fact = round(matrix(c(out.cov$par[c(15,17,17,16)]), ncol=2),3),
            out.std = round(data.frame(loads.cor,disturbs.cor,Rsq=(1-disturbs.cor)),3), cor.fact = round(matrix(c(1,out.cor$par[c(17,17)],1), ncol=2),3),
            fit=data.frame(ll = out.cov$value, 
                     AIC= -2*out.cov$value + 2 * (length(par.init.cov)+ncol(y)), 
                     BIC= -2*out.cov$value + log(nrow(y)) * (length(par.init.cov)+ncol(y)))  #note inclusion of intercepts for total number of par
            )
output

######################
### Confirm lavaan ###
######################

library(lavaan)
y = data.frame(y)
model <- ' F1  =~ X1 + X2 + X3 + X4
           F2  =~ X5 + X6 + X7 + X8 '

fit <- cfa(model, data=y)  #std.lv=T for standardized
summary(fit, fit.measures=TRUE)  #note that lavaan does not count the intercepts among the free params for AIC/BIC, so can take out ncol part above if desired



########################
### Confirm in Mplus ###
########################

library(MplusAutomation)
prepareMplusData(data.frame(y),"C:/Users/mclark19/Desktop/factsim.dat")

# MODEL:
#   F1 BY X1-X4;
#   F2 BY X5-X8;
# 
# OUTPUT:
#   STDYX;
#########################################################################################################
### This more or less comes follows Bollen (1989) for maximum likelihood estimation of a confirmatory ###
### factor analysis. In the following example we will examine a situation where there are two         ###
### underlying (correlated) latent variables for 8 observed responses.  The code as is will only      ### 
### work with this toy data set.  Results are checked against the lavaan package.                     ###
#########################################################################################################

###########################
### Create the data set ###
###########################

library(mvtnorm)
library(psych)

set.seed(123)

# loading matrix
lambda = matrix(c(1,.5,.3,.6,0,0,0,0,
                  0,0,0,0,1,.7,.4,.5),
                nrow=2, byrow=T)

# correlation of factors
phi = matrix(c(1,.25,.25,1), nrow=2, byrow=T)  

# factors and some noise
factors = rmvnorm(1000, mean=rep(0,2), sigma=phi, "chol")
e = rmvnorm(1000, sigma=diag(8))

# observed responses
y = 0 + factors%*%lambda + e

# Examine
#dim(y)
describe(y)
round(cor(y), 3)

#see the factor structure
lazerhawk::corrheat(cor(y))

# example exploratory fa
#fa(y, nfactors=2, rotate="oblimin") 


#########################
### Primary Functions ###
#########################

# measurement model, covariance approach
cfa.cov = function (parms, data) {
  # Arguments- parms: initial values (named); data: raw data
  # Extract paramters by name
  
  require(psych) # for tr
  
  l1 = c(1, parms[grep('l1', names(parms))])      # loadings for factor 1
  l2 = c(1, parms[grep('l2', names(parms))])      # loadings for factor 2
  cov0 = parms[grep('cov', names(parms))]         # factor covariance, variances
  
  # Covariance matrix
  S = cov(data)*((nrow(data)-1)/nrow(data))       # ML covariance div by N rather than N-1, the multiplier adjusts
  
  # loading estimates
  lambda = cbind(c(l1, rep(0,length(l2))),
                 c(rep(0,length(l1)), l2)
                 )
  
  # disturbances
  dist.init = parms[grep('dist', names(parms))]    
  disturbs = diag(dist.init)
  
  # factor correlation
  phi.init = matrix(c(cov0[1], cov0[2], cov0[2], cov0[3]), 2, 2)  #factor cov/correlation matrix
  
  # other calculations and log likelihood
  sigtheta = lambda%*%phi.init%*%t(lambda) + disturbs
  pq = dim(data)[2]  #in Bollen p + q (but for the purposes of this just p) = tr(data)
  #out = -(log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - pq)  #a reduced version; Bollen 1989 p.107
  ll = ((-nrow(data)*pq/2)*log(2*pi)) - (nrow(data)/2)*(log(det(sigtheta)) + tr(S%*%solve(sigtheta)))  #should be same as Mplus H0 loglike
  ll
}

# correlation approach for standardized output; lines correspond to those in cfa.cov 
cfa.cor = function (parms, data) {
  require(psych)
  
  l1 = parms[grep('l1', names(parms))]      # loadings for factor 1
  l2 = parms[grep('l2', names(parms))]      # loadings for factor 2
  cor0 = parms[grep('cor', names(parms))]   # factor correlation
  
  S = cor(data)
  
  lambda = cbind(c(l1, rep(0,length(l2))),
                 c(rep(0,length(l1)), l2)
                 )
  
  dist.init = parms[grep('dist', names(parms))]
  disturbs = diag(dist.init)
  
  phi.init = matrix(c(1, cor0, cor0, 1), ncol=2)
  
  sigtheta = lambda%*%phi.init%*%t(lambda) + disturbs
  pq = dim(data)[2]
  #out = (log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - pq )
  out = ((-nrow(data)*pq/2)*log(2*pi)) - (nrow(data)/2)*(log(det(sigtheta)) + tr(S%*%solve(sigtheta)))
  out
}


####################
### optimization ###
####################
### raw
# initial values
par.init.cov=c(rep(1,6), rep(.05,8), rep(.5,3)) 
names(par.init.cov) = rep(c('l1','l2', 'dist', 'cov'), c(3,3,8,3))

# estimate and extract
out.cov = optim(par=par.init.cov, fn=cfa.cov, data=y, method="L-BFGS-B", lower=0, control=list(fnscale=-1)) 
loads.cov= data.frame(f1=c(1,out.cov$par[1:3], rep(0,4)), f2=c(rep(0,4), 1, out.cov$par[4:6]))
disturbs.cov = out.cov$par[7:14]

# standardized
par.init.cor=c(rep(1,8), rep(.05,8), 0) #for cor
names(par.init.cor) = rep(c('l1','l2', 'dist', 'cor'), c(4,4,8,1))
out.cor = optim(par=par.init.cor, fn=cfa.cor, data=y, method="L-BFGS-B", lower=0, upper=1, control=list(fnscale = -1)) 
loads.cor= matrix(c(out.cor$par[1:4], rep(0,4), rep(0,4), out.cor$par[5:8]), ncol=2)
disturbs.cor = out.cor$par[9:16]


#################################
### Gather output for summary ###
#################################

output=list(raw=list(loadings = round(data.frame(loads.cov, Variances=disturbs.cov), 3), 
            cov.fact = round(matrix(c(out.cov$par[c(15,16,16,17)]), ncol=2), 3)),
            
            standardized=list(loadings = round(data.frame(loads.cor, Variances=disturbs.cor, Rsq=(1-disturbs.cor)), 3), 
            cor.fact = round(matrix(c(1, out.cor$par[c(17,17)], 1), ncol=2), 3)),
            
            fit = data.frame(ll = out.cov$value, 
                                  AIC= -2*out.cov$value + 2 * (length(par.init.cov)+ncol(y)), 
                                  BIC= -2*out.cov$value + log(nrow(y)) * (length(par.init.cov)+ncol(y)))  #note inclusion of intercepts for total number of par
            )
output


###########################
### Confirm with lavaan ###
###########################

library(lavaan)
y = data.frame(y)
model <- ' F1  =~ X1 + X2 + X3 + X4
           F2  =~ X5 + X6 + X7 + X8 '

fit <- cfa(model, data=y, mimic='Mplus', estimator='ML')
fit.std <- cfa(model, data=y, mimic='Mplus', estimator='ML', std.lv=T, std.ov=T) # for standardized

# note that lavaan does not count the intercepts among the free params for AIC/BIC 
# by default, but the mimic='Mplus' should have them correspond to optim's output
summary(fit, fit.measures=TRUE)   


########################
### Confirm in Mplus ###
########################
# If you have access to Mplus you can use Mplus Automation to prepare the data.
# The subsequent code is in Mplus syntax and will produce the above model.
# library(MplusAutomation)
# prepareMplusData(data.frame(y), "factsim.dat")

# MODEL:
#   F1 BY X1-X4;
#   F2 BY X5-X8;
# 
# OUTPUT:
#   STDYX;
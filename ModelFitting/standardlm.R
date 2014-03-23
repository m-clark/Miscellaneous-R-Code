#---------------------------------------------------------------------------------#
# A standard regression model via maximum likelihood or least squares loss.  Also #
# examples for qr decomposition and normal equations. Can serve as an entry point #
# for those starting out to the wider world of computational statistics as        #
# maximum likelihood is the fundamental approach used in most applied statistics, #
# but which is also a key aspect of the Bayesian approach.  Least squares loss    #
# is not confined to the standard lm setting, but is widely used in more          #
# predictive/'algorithmic' approaches e.g. in machine learning and elsewhere.     #
#---------------------------------------------------------------------------------#


##############
# Data Setup #
##############
set.seed(123)  # ensures replication

# predictors and response
N = 100 # sample size
k = 2   # number of desired predictors
X = matrix(rnorm(N*k), ncol=k)  
y = -.5 + .2*X[,1] + .1*X[,2] + rnorm(N, sd=.5)  # increasing N will get estimated values closer to these

dfXy = data.frame(X,y)



#############
# Functions #
#############
# A maximum likelihood approach
lmfuncML = function(par, X, y){
  # arguments- par: parameters to be estimated; X: predictor matrix with intercept 
  # column; y: response
  
  # setup
  beta = par[-1]                               # coefficients
  sigma2 = par[1]                              # error variance
  sigma = sqrt(sigma2)
  N = nrow(X)
  
  # linear predictor
  LP = X%*%beta                                # linear predictor
  mu = LP                                      # identity link in the glm sense
  
  # calculate likelihood
  L = dnorm(y, mean=mu, sd=sigma, log=T)       # log likelihood
#   L =  -.5*N*log(sigma2) - .5*(1/sigma2)*crossprod(y-mu)    # alternate log likelihood form

  -sum(L)                                      # optim by default is minimization, and we want to maximize the likelihood 
                                               # (see also fnscale in optim.control)
}

# An equivalent approach via least squares loss function
lmfuncLS = function(par, X, y){
  # arguments- par: parameters to be estimated; X: predictor matrix with intercept 
  # column, y: response
  
  # setup
  beta = par                                   # coefficients
  
  # linear predictor
  LP = X%*%beta                                # linear predictor
  mu = LP                                      # identity link
  
  # calculate least squares loss function
  L = crossprod(y-mu)
}



##############################
### Obtain Model Estimates ###
##############################
# Setup for use with optim
X = cbind(1, X)

# initial values; note we'd normally want to handle the sigma differently as
# it's bounded by zero, but we'll ignore for demonstration; also sigma2 is not
# required for the LS approach as it is the objective function.
init = c(1, rep(0, ncol(X)));  names(init)=c('sigma2', 'intercept','b1', 'b2')

optlmML = optim(par=init, fn=lmfuncML, X=X, y=y, control=list(reltol=1e-8))
optlmLS = optim(par=init[-1], fn=lmfuncLS, X=X, y=y, control=list(reltol=1e-8))
parsML = optlmML$par
parsLS = c(sigma2=optlmLS$value/(N-k-1), optlmLS$par)  # calculate sigma2 and add



##################
### Comparison ###
##################
### compare to lm which uses QR decomposition
modlm = lm(y~., dfXy)

# Example
# QRX = qr(X)
# Q = qr.Q(QRX)
# R = qr.R(QRX)
# Bhat = solve(R) %*% crossprod(Q, y)
# alternate: qr.coef(QRX, y)

round(rbind(parsML,
            parsLS,
            modlm = c(summary(modlm)$sigma^2, coef(modlm))), 3)

# the slight difference in sigma is roughly maxlike dividing by N vs. N-k-1 in the
# traditional least squares approach; diminishes with increasing N as both tend 
# toward whatever sd^2 you specify when creating the y response above.


### compare to glm; by default assumes gaussian family with identity link and uses lm.fit
modglm = glm(y~., data=dfXy)
summary(modglm)


### via normal equations
coefs = solve(t(X)%*%X) %*% t(X)%*%y  # coefficients

# compare
sqrt(crossprod(y-X%*%coefs)/(N-k-1)); summary(modlm)$sigma; sqrt(modglm$deviance/modglm$df.residual) 
c(sqrt(parsML[1]), sqrt(parsLS[1]))

# rerun by adding 3-4 zeros to the N
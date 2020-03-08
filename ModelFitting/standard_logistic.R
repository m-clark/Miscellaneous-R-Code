#---------------------------------------------------------------------------------#
# A standard logistic regression model via maximum likelihood or exponential
# loss. Can serve as an entry point for those starting out to the wider world of
# computational statistics as maximum likelihood is the fundamental approach used
# in most applied statistics, but which is also a key aspect of the Bayesian
# approach.  Exponential loss is not confined to the standard glm setting, but
# is widely used in more predictive/'algorithmic' approaches e.g. in
# machine learning and elsewhere.

# This follows the standard_lm.R script.
#---------------------------------------------------------------------------------#


##############
# Data Setup #
##############

set.seed(1235)  # ensures replication


# predictors and target

N = 10000 # sample size
k = 2     # number of desired predictors
X = matrix(rnorm(N * k), ncol = k)

# the linear predictor
lp = -.5 + .2 * X[, 1] + .1 * X[, 2] # increasing N will get estimated values closer to these

y = rbinom(N, size = 1, prob = plogis(lp))

dfXy = data.frame(X, y)



#############
# Functions #
#############

# A maximum likelihood approach

logreg_ML = function(par, X, y){
  # arguments- par: parameters to be estimated; X: predictor matrix with intercept 
  # column; y: response
  
  # setup
  beta = par                                # coefficients
  N = nrow(X)
  
  # linear predictor
  LP = X %*% beta                           # linear predictor
  mu = plogis(LP)                           # logit link
  
  # calculate likelihood
  L = dbinom(y, size = 1, prob = mu, log = TRUE)         # log likelihood
  #   L =  y*log(mu) + (1 - y)*log(1-mu)    # alternate log likelihood form
  
  -sum(L)                                   # optim by default is minimization, and we want to maximize the likelihood 
  # (see also fnscale in optim.control)
}

# An equivalent approach via exponential loss function

logreg_exp = function(par, X, y){
  # arguments- par: parameters to be estimated; X: predictor matrix with intercept 
  # column, y: response
  
  # setup
  beta = par                                   # coefficients
  
  # linear predictor
  LP = X %*% beta                              # linear predictor

  # calculate exponential loss function (convert y to -1:1 from 0:1)
  L = sum(exp(-ifelse(y, 1, -1) * .5 * LP))
}


##############################
### Obtain Model Estimates ###
##############################

# Setup for use with optim

X = cbind(1, X)

# initial values

init = rep(0, ncol(X))
names(init)=c('intercept','b1', 'b2')

optlmML = optim(
  par = init,
  fn = logreg_ML,
  X = X,
  y = y,
  control = list(reltol = 1e-8)
)

optglmClass = optim(
  par = init,
  fn = logreg_exp,
  X = X,
  y = y, 
  control = list(reltol = 1e-15)
)

pars_ML  = optlmML$par
pars_exp  = optglmClass$par

##################
### Comparison ###
##################

### compare to glm

modglm = glm(y~., dfXy, family = binomial)

rbind(
  pars_ML,
  pars_exp,
  pars_GLM = coef(modglm)
)

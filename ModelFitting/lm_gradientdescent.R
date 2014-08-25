# Gradient descent for standard linear model.  The function takes arguments
# starting points for the parameters to be estimated, a tolerance to provide a
# stopping point, stepsize, whether to print out iterations, and whether to plot
# the loss over each estimate
# 

##################
### Data Setup ###
##################

# basic data for standard regression. 
set.seed(8675309)
n = 1000
x1 = rnorm(n)
x2 = rnorm(n)
y = 1 + .5*x1 + .2*x2 + rnorm(n)

X = cbind(Intercept=1, x1, x2)


##################################
### gradient descent algorithm ###
##################################

gd = function(par, X, y, tolerance=1e-3, stepsize=1e-3, verbose=T, plotLoss=T){
  # initialize
  beta = par; names(beta) = colnames(X)
  loss = crossprod(X%*%beta - y)
  tol = 1
  iter = 0
  
  while(tol > tolerance){
    LP = X%*%beta
    grad = t(X) %*% (LP - y)
    betaCurrent = beta - stepsize * grad
    tol = max(abs(betaCurrent - beta))
    beta = betaCurrent
    loss = append(loss, crossprod(LP-y))
    iter = iter + 1
    if(verbose && iter%%10 == 0) message(paste('Iteration:', iter))
  }
  
  if(plotLoss) plot(loss, type='l')
  
  list(par=beta, loss=loss, RSE=sqrt(crossprod(LP-y)/(nrow(X)-ncol(X))), 
       iter=iter, fitted=LP)
}

###########
### Run ###
###########

### starting values
init = rep(0, 3)

# for any particular data you'd have to fiddle with the stepsize, which would
# typically be assessed via cross-validation
out = gd(init, X=X, y=y, tolerance = 1e-5, stepsize=.0001)
str(out)

summary(lm(y~x1+x2))


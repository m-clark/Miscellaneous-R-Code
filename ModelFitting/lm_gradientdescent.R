# Gradient descent for standard linear model.  The function takes arguments 
# starting points for the parameters to be estimated, a tolerance or maximum
# iteration value to provide a stopping point, stepsize (or starting stepsize
# for adaptive approach), whether to print out iterations, and whether to plot 
# the loss over each iteration.


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

gd = function(par, X, y, tolerance=1e-3, maxit=1000, stepsize=1e-3, adapt=F,
              verbose=T, plotLoss=T){
  # initialize
  beta = par; names(beta) = colnames(X)
  loss = crossprod(X%*%beta - y)
  tol = 1
  iter = 1
  
  while(tol > tolerance && iter < maxit){
    LP = X%*%beta
    grad = t(X) %*% (LP - y)
    betaCurrent = beta - stepsize * grad
    tol = max(abs(betaCurrent - beta))
    beta = betaCurrent
    loss = append(loss, crossprod(LP-y))
    iter = iter + 1
    if(adapt) stepsize = ifelse(loss[iter] < loss[iter-1],  stepsize*1.2, stepsize*.8)
    if(verbose && iter%%10 == 0) message(paste('Iteration:', iter))
  }
  
  if(plotLoss) plot(loss, type='l', bty='n')
  
  list(par=beta, loss=loss, RSE=sqrt(crossprod(LP-y)/(nrow(X)-ncol(X))), 
       iter=iter, fitted=LP)
}

###########
### Run ###
###########

### starting values
init = rep(0, 3)

# for any particular data you'd have to fiddle with the stepsize, which could 
# be assessed via cross-validation, or alternatively one can use an
# adaptive approach, a simple one of which is implemented in the function
out = gd(init, X=X, y=y, tolerance = 1e-5, stepsize=.0001, adapt=T)
str(out)
round(out$par,5)
summary(lm(y ~ x1 + x2))


# online learning via online gradient descent example
# 

##################
### Data Setup ###
##################

# basic data for standard regression. Note that we are assuming the observations
# come to us over time rather than as a single batch.
set.seed(1234)
x1 = rnorm(1000)
x2 = rnorm(1000)
y = 5 + .5*x1 + .2*x2 + rnorm(1000)

X = cbind(Intercept=1, x1, x2)

#########################################
### online gradient descent algorithm ###
#########################################

ogd = function(starts, X, y, stepsizeTau=0, stepsize=1){
  # initialize
  beta = starts; names(beta) = colnames(X)
  betavec = matrix(0, nrow(X), ncol=length(beta))
  regret = numeric(nrow(X))
  s = 0
  
  for (t in 1:nrow(X)){
    Xt = X[t,, drop=F]
    yt = y[t]
    LP = Xt%*%beta
    grad = stepsize * t(Xt) %*% (LP - yt)
    s = s + grad^2
    beta = beta - stepsize*grad/(stepsizeTau + sqrt(s))  # adagrad
    betavec[t,] = beta
    regret[t] = crossprod(X[1:t,]%*%beta -y[t])/t - ifelse(t=1, 0, min(regret)/t)
  }
  
  LP = X%*%beta
  loss = crossprod(LP - y)
  list(par=beta, parvec=betavec, loss=loss, regret=regret,
       RSE=loss[length(loss)]/nrow(X), fitted=LP)
}


###########
### Run ###
###########

### starting values
init = rep(0, 3)
# debugonce(ogd)
# for any particular data you'd have to fiddle with the stepsize, which would
# typically be assessed via cross-validation
out = ogd(init, X=X, y=y, stepsizeTau=0, stepsize=1)
str(out)
summary(lm(y~x1+x2))
# plot(out$fitted, y)
plot(out$regret, type='l')

library(ggplot2); library(reshape2)
gd = melt(out$parvec); colnames(gd) = c('T', 'Parameter', 'Value')
ggplot(aes(x=T, y=Value), data=gd) +
  geom_path(aes(group=Parameter)) +
  theme_minimal()

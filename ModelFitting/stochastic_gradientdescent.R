# 'Online' learning via online/stochastic gradient descent.  See also, standard 
# gradient descent in lm_gradientdescent.R In the following, we have basic data 
# for standard regression, but in this 'online' learning case, we are assuming 
# each observation comes to us as a stream over time rather than as a single 
# batch, and would continue coming in.  Note that there are plenty of variations
# of this, an it can be applied in the batch case as well.  Currently no 
# stopping point is implemented in order to trace results over all data
# points/iterations.

##################
### Data Setup ###
##################
set.seed(1234)
n = 1000
x1 = rnorm(n)
x2 = rnorm(n)
y = 1 + .5*x1 + .2*x2 + rnorm(n)

X = cbind(Intercept=1, x1, x2)



#############################################
### stochastic gradient descent algorithm ###
#############################################

sgd = function(par, X, y, stepsizeTau=0, stepsize=1, average=F){
  # initialize
  beta = par; names(beta) = colnames(X)
  betamat = matrix(0, nrow(X), ncol=length(beta))        # Collect all estimates
  fits = NA                                              # fitted values
  s = 0                                                  # for adagrad
  loss = NA                                              # Collect loss at each point
  
  for (t in 1:nrow(X)){
    Xt = X[t,, drop=F]
    yt = y[t]
    LP = Xt%*%beta                                       # matrix operations not necessary, 
    grad =  t(Xt) %*% (LP - yt)                          # but makes consistent standard gd R file
    s = s + grad^2
    beta = beta - stepsize*grad/(stepsizeTau + sqrt(s))  # adagrad approach
    if(average & t>1){
      beta =  beta - 1/t*(betamat[t-1,]-beta) 
    } 
    betamat[t,] = beta
    fits[t] = LP
    loss[t] = (LP-yt)^2
  }
  
  LP = X%*%beta
  lastloss = crossprod(LP - y)
  list(par=beta, parvec=betamat, loss=loss, 
       RSE=sqrt(sum(loss)/(nrow(X)-ncol(X))), fitted=fits)
}



###########
### Run ###
###########

### starting values
init = rep(0, 3)

# for any particular data you might have to fiddle with the stepsize, perhaps
# choosing one based on cross-validation with old data
out = sgd(init, X=X, y=y, stepsizeTau=1e-3, stepsize=.2, average = T)
str(out)
out$par

# compare
summary(lm(y ~ x1 + x2))
coef1 = coef(lm(y ~ x1 + x2))


# visualize estimates
library(ggplot2); library(reshape2)
gd = melt(out$parvec); colnames(gd) = c('T', 'Parameter', 'Value')
gd$Parameter = factor(gd$Parameter, labels=colnames(X))
ggplot(aes(x=T, y=Value, group=Parameter, color=Parameter), data=gd) +
  geom_path(aes()) +
  geom_point(data=gd[gd$T==n,], size=3) +
  geom_text(aes(label=round(Value,2)), hjust=-.5, angle=45, 
            data=gd[gd$T==n,], size=4, show_guide=F) +
  theme_minimal()



######################################
### add alternately generated data ###
######################################
# data shift
set.seed(1234)
n2 = 1000
x1.2 = rnorm(n2)
x2.2 = rnorm(n2)
y2 = -1 + .25*x1.2 - .25*x2.2 + rnorm(n2)

X2 = rbind(X, cbind(1, x1.2, x2.2))
coef2 = coef(lm(y2 ~ x1.2 + x2.2))

y2 = c(y, y2)

n3 = 1000
x1.3 = rnorm(n3)
x2.3 = rnorm(n3)
y3 = 1 - .25*x1.3 + .25*x2.3 + rnorm(n3)
coef3 = coef(lm(y3 ~ x1.3 + x2.3))

X3 = rbind(X2, cbind(1, x1.3, x2.3))
y3 = c(y2, y3)


###########
### Run ###
###########
out2 = sgd(init, X=X3, y=y3, stepsizeTau=1e-3, stepsize=1, average = T)
str(out2)

# compare with lm for each data part
out2$parvec[c(n, n+n2, n+n2+n3),]
rbind(coef1, coef2, coef3)

# visualize estimates
gd = melt(out2$parvec); colnames(gd) = c('T', 'Parameter', 'Value')
gd$Parameter = factor(gd$Parameter, labels=colnames(X))
ggplot(aes(x=T, y=Value, group=Parameter, color=Parameter), data=gd) +
  geom_path(aes()) +
  geom_point(data=gd[gd$T==(c(n, n+n2, n+n2+n3)),], size=3) +
  geom_text(aes(label=round(Value,2)), hjust=-.5, angle=45, data=gd[gd$T==(c(n, n+n2, n+n2+n3)),], 
            size=4, show_guide=F) +
  theme_minimal()
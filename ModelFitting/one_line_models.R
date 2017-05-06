
# One line models ---------------------------------------------------------

# A terrible idea put to good use in demonstrating vectorization and other 
# efficient coding practices, after data set up, parameters are learned via a 
# one-line minimizing/maximizing function, typically using optim.  Or maybe it's
# just a fun exercise. Who can say?

library(tidyverse)

# Standard Regression -----------------------------------------------------

# data setup
set.seed(8675309)
N = 100
npreds = 5
X = cbind(1, matrix(rnorm(N*npreds), ncol=npreds))
beta = runif(ncol(X), -1, 1)
y = X %*% beta + rnorm(nrow(X))

# Normal equations
crossprod(solve(crossprod(X)), crossprod(X, y))

# use lm for comparison
coef(lm.fit(X,y))

# run model via maximum likelihoood
optim(
  rep(0, ncol(X)),
  fn = function(b, X, y) crossprod(y - X %*% b),  # model function
  X = X,
  y = y,
  method = 'BFGS'
)$par

# use lm for comparison
coef(lm.fit(X,y))



# Logistic Regression -----------------------------------------------------

# data setup
y01 = rbinom(N, size=1, p=plogis(X %*% beta))

# run model
optim(
  rep(0, ncol(X)),
  fn = function(b, X, y) -sum(dbinom(y, size = 1, prob = plogis(X %*% b), log = T)),   # model function
  X = X,
  y = y01,
  method = 'BFGS'
)$par

# use glm for comparison
coef(glm.fit(X, y01, family=binomial()))



# Poisson Regression -----------------------------------------------------

# data setup
y_count = rpois(N, exp(X %*% beta))

# run model
optim(
  rep(0, ncol(X)),
  fn = function(b, X, y) -sum(dpois(y, exp(X %*% b), log = T)),   # model function
  X = X,
  y = y_count,
  method = 'BFGS'
)$par

# use glm for comparison
coef(glm.fit(X, y_count, family=poisson()))



#  Naive bayes for binary data --------------------------------------------

# data setup
x = matrix(sample(0:1, 50, replace = T), ncol=5)
xf = data.frame(lapply(data.frame(x), factor))
y = sample(0:1, 10, prob=c(.25, .75), replace=T)


# use e1071 package for comparison
library(e1071)
m = naiveBayes(xf, y)

# run model
lapply(xf, function(var) t(prop.table(table(' '=var, y), margin=2)))

m



# Cox Regression ----------------------------------------------------------

# data setup
dur = 1:10
kittyblarg= rnorm(10)                                 # something happened to kitty!
kittyhappy = rep(0:1,times=5)                         # is kitty happy?
kittydied = sample(0:1,10,replace=T)                  # kitty died! oh noes!
d = data.frame(kittyblarg,kittyhappy,dur,kittydied)[order(dur),]

X = cbind(kittyblarg, kittyhappy)

# run model
optim(par = rep(0, ncol(X)), 
      fn = function(b, X, died, t) -sum(sapply(1:nrow(X), function(i) died[i]*((X%*%b)[i] - log(sum(exp((X%*%b)[1:nrow(X)>=i])))))), 
      X = X,
      died = d$kittydied, 
      t=dur, 
      method="BFGS")$par

# use survival package for comparison
library(survival)
coef(coxph(Surv(dur, kittydied) ~ kittyblarg + kittyhappy))




# Mixed Model -------------------------------------------------------------

# data setup
data(sleepstudy, package='lme4')
X = model.matrix(~Days, sleepstudy)
Z = model.matrix(~factor(sleepstudy$Subject)-1)
colnames(Z) =  paste0('Subject_', unique(sleepstudy$Subject))  # for cleaner presentation later
rownames(Z) = paste0('Subject_', sleepstudy$Subject)
y = sleepstudy$Reaction


# run model
optim(par = c(0, 0), 
      fn = function(b, X, Z, y) -mvtnorm::dmvnorm(y, X%*%coef(lm.fit(X, y)), tcrossprod(Z)*exp(b[1])^2 + diag(nrow(X))*exp(b[2])^2, log=T),
      X = X,
      Z = Z,
      y = y,
      method='L-BFGS',
      lower=c(0,0))$par %>% exp

# check with lme4
library(lme4)
lmemod = lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
VarCorr(lmemod)




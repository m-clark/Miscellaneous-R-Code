# A very simple implementation of extreme learning machine for regression for 
# demonstration. See elmNN and ELMR for some R package implementations. I add 
# comparison to generalized additive models (elm/neural networks and GAMs are
# adaptive basis function models).

# http://www.extreme-learning-machines.org
# G.-B. Huang, Q.-Y. Zhu and C.-K. Siew, "Extreme Learning Machine: Theory and Applications"

elm <- function(X, y, n_hidden=NULL, active_fun=tanh) {
  # X: an N observations x p features matrix
  # y: the target
  # n_hidden: the number of hidden nodes
  # active_fun: activation function
  pp1 = ncol(X) + 1
  w0 = matrix(rnorm(pp1*n_hidden), pp1, n_hidden)       # random weights
  h = active_fun(cbind(1, scale(X)) %*% w0)             # compute hidden layer
  B = MASS::ginv(h) %*% y                               # find weights for hidden layer
  fit = h %*% B                                         # fitted values
  list(fit= fit, loss=crossprod(fit - y), B=B, w0=w0)
}



# one variable, complex function -------------------------------------------
library(tidyverse); library(mgcv)
set.seed(123)
n = 5000
x = runif(n)
# x = rnorm(n)
mu = sin(2*(4*x-2)) + 2*exp(-(16^2)*((x-.5)^2))
y = rnorm(n, mu, .3)
# qplot(x, y)
d = data.frame(x,y) 

X_ = as.matrix(x, ncol=1)

test = elm(X_, y, n_hidden=100)
str(test)
# qplot(x, y) + geom_line(aes(y=test$fit), color='#1e90ff')
cor(test$fit[,1], y)^2

gam_comparison = gam(y~s(x))
summary(gam_comparison)$r.sq


d %>% 
  mutate(fit_elm = test$fit,
         fit_gam = fitted(gam_comparison)) %>% 
  ggplot() + 
  geom_point(aes(x, y), alpha=.1) +
  geom_line(aes(x, y=fit_elm), color='#1e90ff') + 
  geom_line(aes(x, y=fit_gam), color='darkred')




# motorcycle accident data ------------------------------------------------

data('mcycle', package='MASS')
x = mcycle[,1]
X_ = matrix(x, ncol=1)
y = mcycle[,2]

test = elm(X_, y, n_hidden=100)
cor(test$fit[,1], y)^2

gam_comparison = gam(y~s(x))
summary(gam_comparison)$r.sq

qplot(x, y) +
  geom_line(aes(y=test$fit), color='#1e90ff') + 
  geom_line(aes(y=fitted(gam_comparison)), color='darkred')




# add covariates ----------------------------------------------------------

d = gamSim(eg=7, n=10000)
X_ = as.matrix(d[,2:5])
y = d[,1]

n_nodes = c(10, 25, 100, 250, 500, 1000)
test = lapply(n_nodes, function(n) elm(X_, y, n_hidden=n))       # this will take a few seconds
final_n = which.min(sapply(test, function(x) x$loss))
best = test[[final_n]]
# str(best)
qplot(best$fit[,1], y, alpha=.2)
cor(best$fit[,1], y)^2

gam_comparison = gam(y~s(x0) + s(x1) + s(x2) + s(x3), data=d)
gam.check(gam_comparison)
summary(gam_comparison)$r.sq


test_data0 = gamSim(eg=7)  # default n = 400
test_data =  cbind(1, scale(test_data0[,2:5]))

elm_prediction = tanh(test_data %*% best$w0) %*% best$B          # remember to use your specific activation function here
gam_prediction = predict(gam_comparison, newdata=test_data0)
cor(data.frame(elm_prediction, gam_prediction), test_data0$y)^2

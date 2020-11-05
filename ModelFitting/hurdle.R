#' ---
#' title: "Hurdle Models"
#' author: "Michael Clark"
#' date: ""
#' ---
#' 
#' 
#' # Poisson
hurdpoisloglik = function(y, X, par) {
  # Extract parameters
  logitpars = par[grep('logit', names(par))]
  poispars  = par[grep('pois', names(par))]
  
  # Logit model part
  Xlogit = X
  ylogit = ifelse(y == 0, 0, 1)
  
  LPlogit = Xlogit %*% logitpars
  mulogit = plogis(LPlogit)
  
  # Calculate the likelihood
  logliklogit = -sum( ylogit*log(mulogit) + (1 - ylogit)*log(1 - mulogit) )  
  
  # Poisson part
  Xpois = X[y > 0, ]
  ypois = y[y > 0]
  
  mupois = exp(Xpois %*% poispars)
  
  # Calculate the likelihood
  loglik0    = -mupois
  loglikpois = -sum(dpois(ypois, lambda = mupois, log = TRUE)) + sum(log(1 - exp(loglik0)))
  
  # combine likelihoods
  loglik = loglikpois + logliklogit
  loglik
}


hurdNBloglik = function(y, X, par) {
  # Extract parameters
  logitpars  = par[grep('logit', names(par))]
  NegBinpars = par[grep('NegBin', names(par))]
  
  theta = exp(par[grep('theta', names(par))])
  
  # Logit model part
  Xlogit = X
  ylogit = ifelse(y == 0, 0, 1)
  
  LPlogit = Xlogit%*%logitpars
  mulogit =  plogis(LPlogit)
  
  # Calculate the likelihood
  logliklogit = -sum( ylogit*log(mulogit) + (1 - ylogit)*log(1 - mulogit) )
  
  #NB part
  XNB = X[y > 0, ]
  yNB = y[y > 0]
  
  muNB = exp(XNB %*% NegBinpars)
  
  # Calculate the likelihood
  loglik0  = dnbinom(0,   mu = muNB, size = theta, log = TRUE)
  loglik1  = dnbinom(yNB, mu = muNB, size = theta, log = TRUE)
  loglikNB = -( sum(loglik1) - sum(log(1 - exp(loglik0))) )
  
  # combine likelihoods
  loglik = loglikNB + logliklogit
  loglik
}


#' # Data Import

#' Import a simple data set. Example from the Stata help file for zinb command;
#' can compare results with hnblogit command

library(haven)

fish = read_dta("http://www.stata-press.com/data/r11/fish.dta")

# Get some starting values.

init_mod = glm(
  count ~ persons + livebait,
  data = fish,
  family = poisson,
  x = TRUE,
  y = TRUE
)

#' for these functions, a named vector for the starting values
starts = c(logit = coef(init_mod), pois = coef(init_mod))  

#' # Poisson hurdle 

#' Use `optim` to estimate parameters. I fiddle with some options to reproduce the
#' hurdle function as much as possible.
#' 
optPois1 = optim(
  par = starts,
  fn  = hurdpoisloglik,
  X   = init_mod$x,
  y   = init_mod$y,
  control = list(maxit = 5000, reltol = 1e-12),
  hessian = TRUE
)
# optPois1

#' Extract the elements from the output to create a summary table.
B  = optPois1$par
se = sqrt(diag(solve(optPois1$hessian)))
Z  = B/se
p  = ifelse(Z >= 0, pnorm(Z, lower = FALSE)*2, pnorm(Z)*2)
summarytable = round(data.frame(B, se, Z, p), 3)

list(summary = summarytable, ll = optPois1$value)

#' Compare to hurdle from pscl package.
library(pscl)

poismod = hurdle(
  count ~ persons + livebait,
  data = fish,
  zero.dist = "binomial",
  dist = "poisson"
)

summary(poismod)


#' # Negative Binomial hurdle


starts =  c(
  logit  = coef(init_mod),
  NegBin = coef(init_mod),
  theta  = 1
)

optNB1 = optim(
  par = starts,
  fn  = hurdNBloglik,
  X   = init_mod$x,
  y   = init_mod$y,
  control = list(maxit = 5000, reltol = 1e-12),
  method  = "BFGS",
  hessian = TRUE
)
# optNB1 

B  = optNB1$par
se = sqrt(diag(solve(optNB1$hessian)))
Z  = B/se
p  = ifelse(Z >= 0, pnorm(Z, lower = FALSE)*2, pnorm(Z)*2)
summarytable = round(data.frame(B, se, Z, p), 3)
list(summary = summarytable, ll = optNB1$value)

NBmod = hurdle(
  count ~ persons + livebait,
  data = fish,
  zero.dist = "binomial",
  dist = "negbin"
)

summary(NBmod)


#' ---
#' title: "Zero-inflated Negative Binomial Model"
#' author: "Michael Clark"
#' date: ""
#' ---
#' 
#' 
#' Log likelihood function to estimate parameters for a Zero-inflated Negative Binomial model. With
#' examples and comparison to pscl package output. Also includes approach based on Hilbe GLM text.
#' see also: https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/poiszeroinfl.R

ZINB = function(y, X, par) {
  # arguments are response y, predictor matrix X, and parameter named starting points of 'logit', 'negbin', and 'theta'
  
  # Extract parameters
  logitpars  = par[grep('logit', names(par))]
  negbinpars = par[grep('negbin', names(par))]
  theta = exp(par[grep('theta', names(par))])
  
  # Logit part; in this function Xlogit = Xnegbin but one could split X argument into Xlogit and Xnegbin for example
  Xlogit  = X
  LPlogit = Xlogit %*% logitpars
  logi0   = plogis(LPlogit) 
  
  # Negbin part
  Xnegbin = X
  munb = exp(Xnegbin %*% negbinpars)
  
  # LLs
  logliklogit  = log( logi0 + exp(log(1 - logi0) + suppressWarnings(dnbinom(0, size = theta, mu = munb, log = TRUE))) )
  logliknegbin = log(1 - logi0) + suppressWarnings(dnbinom(y, size = theta, mu = munb, log = TRUE))
  
  # Hilbe formulation
  # theta part 
  # alpha = 1/theta  
  # m = 1/alpha
  # p = 1/(1 + alpha*munb)
  
  # logliklogit = log( logi0 + (1 - logi0)*(p^m) )
  # logliknegbin = log(1-logi0) + log(gamma(m+y)) - log(gamma(m)) + m*log(p) + y*log(1-p)   # gamma(y+1) not needed
  
  y0 = y == 0   # 0 values
  yc = y > 0    # Count part
  
  loglik = sum(logliklogit[y0]) + sum(logliknegbin[yc])
  -loglik
}


#' Get the data
library(haven)

fish = read_dta("http://www.stata-press.com/data/r11/fish.dta")


#' Get starting values or simply do zeros
#' for this function, a named vector for the starting values
#' for zinb: 'logit', 'negbin', 'theta'
init.mod = model.matrix(count ~ persons + livebait, data = fish) # to get X matrix

startlogi  = glm(count == 0 ~ persons + livebait, data = fish, family = "binomial")
startcount = glm(count ~ persons + livebait, data = fish, family = "poisson")

starts = c(
  negbin = coef(startcount),
  logit = coef(startlogi),
  theta = 1
)  
# starts = c(negbin = rep(0, 3),
#            logit = rep(0, 3),
#            theta = log(1))


#' Estimate with optim function
optNB1 = optim(
  par = starts ,
  fn  = ZINB,
  X   = init.mod,
  y   = fish$count,
  method  = "BFGS",
  control = list(maxit = 5000, reltol = 1e-12),
  hessian = TRUE
)
# optNB1

#' Comparison
# Extract for clean display
B  = optNB1$par
se = sqrt(diag(solve((optNB1$hessian))))
Z  = B/se
p  = pnorm(abs(Z), lower = FALSE)*2

# pscl results
library(pscl)
zinbmod1 = zeroinfl(count ~ persons + livebait, data = fish, dist = "negbin")
summary(zinbmod1)
round(data.frame(B, se, Z, p), 4)  # note that theta here is actually log(theta)


#' Optional data set


#' Get the data
data("bioChemists", package = "pscl")

#' Get starting values or simply do zeros
init.mod   = model.matrix(art ~ fem + mar + kid5 + phd + ment, data = bioChemists) # to get X matrix
startlogi  = glm(art==0 ~ fem + mar + kid5 + phd + ment, data = bioChemists, family = "binomial")
startcount = glm(art ~ fem + mar + kid5 + phd + ment, data = bioChemists, family = "quasipoisson")

starts = c(
  negbin = coef(startcount),
  logit  = coef(startlogi),
  theta  = summary(startcount)$dispersion
)  
# starts = c(negbin = rep(0, 6),
#            logit = rep(0, 6),
#            theta = 1)


#' Estimate with optim function
optNB2 = optim(
  par = starts ,
  fn  = ZINB,
  X   = init.mod,
  y   = bioChemists$art,
  method  = "BFGS",
  control = list(maxit = 5000, reltol = 1e-12),
  hessian = TRUE
)
# optNB2


#' Comparison
# Extract for clean display.
B  = optNB2$par
se = sqrt(diag(solve((optNB2$hessian))))
Z  = B/se
p  = pnorm(abs(Z), lower = FALSE)*2

# pscl results
library(pscl)
zinbmod = zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")
summary(zinbmod)
round(data.frame(B,se, Z, p), 4)
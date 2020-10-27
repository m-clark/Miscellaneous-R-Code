########################################################################################################
### Log likelihood function to estimate parameters for a Zero-inflated Poisson model. With examples  ###
### and comparison to pscl package output. Also includes approach based on Hilbe GLM text.           ###
### see also: https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/NBzeroinfl.R ###
########################################################################################################

ZIP = function(y, X, par) {
  # arguments are response y, predictor matrix X, and parameter named starting points of 'logit' and 'pois'
  
  # Extract parameters
  logitpars = par[grep('logit', names(par))]   
  poispars  = par[grep('pois', names(par))]     
  
  # Logit part; in this function Xlogit = Xpois but one could split X argument into Xlogi and Xpois for example
  Xlogit  = X
  LPlogit = Xlogit %*% logitpars
  logi0   = plogis(LPlogit)  # alternative 1/(1+exp(-LPlogit))
    
  # Poisson part
  Xpois  = X
  mupois = exp(Xpois %*% poispars)
  
  # LLs
  logliklogit = log( logi0 + exp(log(1 - logi0) - mupois) )
  loglikpois  = log(1 - logi0) + dpois(y, lambda = mupois, log = TRUE)
  
  # Hilbe formulation
  # logliklogit = log(logi0 + (1 - logi0)*exp(- mupois) )
  # loglikpois = log(1-logi0) -mupois + log(mupois)*y     #not necessary: - log(gamma(y+1))
    
  y0 = y == 0  # 0 values
  yc = y > 0   # Count part

  loglik = sum(logliklogit[y0]) + sum(loglikpois[yc])
  -loglik
}

### Get the data ###
library(haven)
library(pscl)

fish = read_dta("http://www.stata-press.com/data/r11/fish.dta")


### Get starting values or simply do zeros ###
# for this function, a named vector for the starting values; for zip: need 'logit', 'pois'
init.mod = glm(
  count ~ persons + livebait,
  data = fish,
  x = TRUE,
  y = TRUE,
  "poisson"
)

# starts = c(logit=coef(init.mod), pois=coef(init.mod))  
starts = c(rep(0, 3), rep(0, 3))
names(starts) = c(paste0('pois.', names(coef(init.mod))), paste0('logit.', names(coef(init.mod))))
                                                               

### Estimate with optim function ###
optPois1 = optim(
  par = starts ,
  fn  = ZIP,
  X   = init.mod$x,
  y   = init.mod$y,
  method  = "BFGS",
  control = list(maxit = 5000, reltol = 1e-12),
  hessian = TRUE
)

# optPois1


### Comparison ###
# Extract for clean display
B = optPois1$par
se = sqrt(diag(solve((optPois1$hessian))))
Z = B/se
p = pnorm(abs(Z), lower = FALSE)*2

# pscl results
zipoismod = zeroinfl(count ~ persons + livebait, data = fish, dist = "poisson") 
summary(zipoismod)
round(data.frame(B, se, Z, p), 4)

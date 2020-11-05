#' ---
#' title: " GLM estimation"
#' subtitle: "Newton and IRLS"
#' author: "Michael Clark"
#' css: '../other.css'
#' highlight: pygments
#' date: ""
#' ---
#'
#' # GLM estimation examples

#' Examples of maximum likelihood estimation via a variety of means.  See the
#' gradientdescent.R script for that approach.  Here we demonstrate Newton's and
#' Iterated Reweighted Least Squares approaches via logistic regression.
#' 
#' 
#' For the following, I had Murphy's PML text open and more or less followed the 
#' algorithms in chapter 8.  Note that for Newton's method, this doesn't
#' implement a line search to find a more optimal stepsize at a given iteration.
#' 
#' # Data Prep
#' 
#' Predict graduate school admission based on gre, gpa, and school rank
#' (higher=more prestige). See corresponding demo here:
#' https://stats.idre.ucla.edu/stata/dae/logistic-regression/. The only
#' difference is that I treat rank as numeric rather than categorical.


admit = haven::read_dta('https://stats.idre.ucla.edu/stat/stata/dae/binary.dta')

comparison_model = glm(admit ~ gre + gpa + rank, data = admit, family = binomial)

summary(comparison_model)

X = model.matrix(comparison_model)
y = comparison_model$y


#' # Newton's method

newton <- function(
  X,
  y,
  tol  = 1e-12,
  iter = 500,
  stepsize = .5
  ) {
  
  # Args: 
  # X: model matrix
  # y: target
  # tol: tolerance
  # iter: maximum number of iterations
  # stepsize: (0, 1)
  
  # intialize
  int     = log(mean(y) / (1 - mean(y)))         # intercept
  beta    = c(int, rep(0, ncol(X) - 1))
  currtol = 1
  it = 0
  ll = 0
  
  while (currtol > tol && it < iter) {
    it = it +1
    ll_old = ll
    
    mu = plogis(X %*% beta)[,1]
    g  = crossprod(X, mu-y)               # gradient
    S  = diag(mu*(1-mu)) 
    H  = t(X) %*% S %*% X                 # hessian
    beta = beta - stepsize * solve(H) %*% g

    ll = sum(dbinom(y, prob = mu, size = 1, log = TRUE))
    currtol = abs(ll - ll_old)
  }
  
  list(
    beta = beta,
    iter = it,
    tol  = currtol,
    loglik = ll
  )
}


newton_result = newton(
  X = X,
  y = y,
  stepsize = .9,
  tol = 1e-8      # tol set to 1e-8 as in glm default
) 

newton_result
comparison_model

rbind(
  newton = unlist(newton_result),
  glm_default = c(
    beta = coef(comparison_model),
    comparison_model$iter,
    tol = NA,
    loglik = -logLik(comparison_model)
  )
)


#' # IRLS 
#' Note that `glm` is actually using IRLS, so the results from this should be
#' fairly spot on.

irls <- function(X, y, tol = 1e-12, iter = 500) {
  
  # intialize
  int  = log(mean(y) / (1 - mean(y)))   # intercept
  beta = c(int, rep(0, ncol(X) - 1))
  currtol = 1
  it = 0
  ll = 0
  
  while (currtol > tol && it < iter) {
    it = it + 1
    ll_old = ll
    
    eta  = X %*% beta
    mu   = plogis(eta)[,1]
    s    = mu * (1 - mu)
    S    = diag(s)
    z    = eta + (y-mu)/s
    beta = solve(t(X) %*% S %*% X) %*% (t(X) %*% (S %*% z))
    
    ll = sum(
      dbinom(
        y,
        prob = plogis(X %*% beta),
        size = 1,
        log = T
      )
    )
    
    currtol = abs(ll - ll_old)
  }
  
  list(
    beta = beta,
    iter = it,
    tol  = currtol,
    loglik  = ll,
    weights = plogis(X %*% beta) * (1 - plogis(X %*% beta))
  )
}

#' `tol` set to 1e-8 as in `glm` default.
irls_result = irls(X = X, y = y, tol = 1e-8) 

str(irls_result)
comparison_model

#' # Comparison
#' 
#' Compare all results.
rbind(
  newton = unlist(newton_result),
  irls   = unlist(irls_result[-length(irls_result)]),
  glm_default = c(
    beta = coef(comparison_model),
    comparison_model$iter,
    tol = NA,
    loglik = logLik(comparison_model)
  )
)


#' compare weights 
head(cbind(irls_result$weights,
           comparison_model$weights))



#' # Source
#' Base R source code found at https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/newton_irls.R
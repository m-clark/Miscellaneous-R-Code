#' ---
#' title: " L1 (lasso) regularization"
#' author: "Michael Clark"
#' css: '../other.css'
#' highlight: pygments
#' date: ""
#' ---
#' 
#' See Tibshirani (1996) for the source, or Murphy PML (2012) for a nice
#' overview (watch for typos in depictions). A more conceptual depiction of the
#' lasso can be found in penalized_ML.R.
#' 
#' # Coordinate descent
#' 
#' 
lasso <- function(
  X,                   # model matrix
  y,                   # target
  lambda  = .1,        # penalty parameter
  soft    = TRUE,      # soft vs. hard thresholding
  tol     = 1e-6,      # tolerance
  iter    = 100,       # number of max iterations
  verbose = TRUE       # print out iteration number
) {
  
  # soft thresholding function
  soft_thresh <- function(a, b) {
    out = rep(0, length(a))
    out[a >  b] = a[a > b] - b
    out[a < -b] = a[a < -b] + b
    out
  }
  
  w = solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,y)
  tol_curr = 1
  J = ncol(X)
  a = rep(0, J)
  c_ = rep(0, J)
  i = 1
  
  while (tol < tol_curr && i < iter) {
    w_old = w 
    a = colSums(X^2)
    l = length(y)*lambda  # for consistency with glmnet approach
    c_ = sapply(1:J, function(j)  sum( X[,j] * (y - X[,-j] %*% w_old[-j]) ))
    if (soft) {
      for (j in 1:J) {
        w[j] = soft_thresh(c_[j]/a[j], l/a[j])
      }
    }
    else {
      w = w_old
      w[c_< l & c_ > -l] = 0
    }
    
    tol_curr = crossprod(w - w_old)  
    i = i + 1
    if (verbose && i%%10 == 0) message(i)
  }
  
  w
}

#' # Data setup
#' 
#' 
set.seed(8675309)
N = 500
p = 10
X = scale(matrix(rnorm(N*p), ncol=p))
b = c(.5, -.5, .25, -.25, .125, -.125, rep(0, p-6))
y = scale(X %*% b + rnorm(N, sd=.5))
lambda = .1


# debugonce(lasso)

#' Note, if `lambda=0`, result is the same as  `lm.fit`.
#' 
#' 
result_soft = lasso(
  X,
  y,
  lambda = lambda,
  tol = 1e-12,
  soft = TRUE
)

result_hard = lasso(
  X,
  y,
  lambda = lambda,
  tol    = 1e-12,
  soft   = FALSE
)




#' `glmnet` is by default a mixture of ridge and lasso penalties, setting alpha
#' = 1 reduces to lasso (alpha=0 would be ridge). We set the lambda to a couple
#' values while only wanting the one set to the same lambda value as above (s).


library(glmnet)

glmnet_res = coef(
  glmnet(
    X,
    y,
    alpha  = 1,
    lambda = c(10, 1, lambda),
    thresh = 1e-12,
    intercept = FALSE
  ),
  s = lambda
)

library(lassoshooting)

ls_res = lassoshooting(
  X = X,
  y = y,
  lambda = length(y) * lambda,
  thr = 1e-12
)


#' # Comparison

data.frame(
  lm = coef(lm(y ~ . - 1, data.frame(X))),
  lasso_soft = result_soft,
  lasso_hard = result_hard,
  lspack = ls_res$coef,
  glmnet = glmnet_res[-1, 1],
  truth  = b
)



#' # Source
#' Base R source code found at https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/lasso.R


#' for some more detailed R code, check out
#' http://jocelynchi.com/a-coordinate-descent-algorithm-for-the-lasso-problem (now defunct, but might find relevant article at the website)


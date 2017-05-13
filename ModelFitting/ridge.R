# l2 (ridge) regularization -----------------------------------------------

# Compare to lasso.R

# 
ridge <- function(w, X, y, lambda = .1) {
  # X: model matrix; y: target; lambda: penalty parameter; w: the weights/coefficients
  
  crossprod(y - X%*%w) + lambda*length(y)*crossprod(w)
}


set.seed(8675309)
N = 500
p = 10
X = scale(matrix(rnorm(N*p), ncol=p))
b = c(.5, -.5, .25, -.25, .125, -.125, rep(0, 4))
y = scale(X %*% b + rnorm(N, sd=.5))


# note, if lambda=0, result is lm.fit
result_ridge = optim(rep(0, ncol(X)), ridge, X=X, y=y, lambda=.1, method='BFGS')
result_ridge2 =  solve(crossprod(X) + diag(length(y)*.1, ncol(X))) %*% crossprod(X,y)



library(glmnet)
# glmnet is by default a mixture of ridge and lasso penalties, setting alpha = 0
# reduces to ridge (alpha=1 would be lasso)
glmnet_res = coef(glmnet(X, y, alpha=0, lambda=c(10, 1, .1), thresh=1e-12, intercept=F), s=.1)



# comparison
data.frame(lm=coef(lm(y~.-1, data.frame(X))),
           ridge=result_ridge$par,
           ridge2=result_ridge2,
           glmnet = glmnet_res[-1, 1],
           truth=b)

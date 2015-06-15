set.seed(3)
n = 20
x1 = runif(n)
x2 = runif(n)
X = matrix(c(x1, x2), ncol = 2) # design matrix
y = 2 + 3 * x1 + rnorm(n, sd = 0.25)
my.inv = function(X, eps = 1e-12) {
  eig.X = eigen(X, symmetric = T)
  P = eig.X[[2]] 
  lambda = eig.X[[1]] 
  ind = lambda > eps
  lambda[ind] = 1/lambda[ind] 
  lambda[!ind] = 0
  ans = P %*% diag(lambda) %*% t(P)
  return(ans)
}
rk = function(s, t) {
  crossprod(s, t)              # MC note: slightly more succinct than the original double loop; also sum(s*t)
} 
get.gramm = function(X, rkfunc=rk) {
  apply(X, 1, function(Row) apply(X, 1, function(tRow) rkfunc(Row, tRow)))  
}
ridge.regression = function(X, y, lambda) {
  Gramm = get.gramm(X)                               # Gramm matrix (nxn)
  n = length(y)                                      # n=length of y
  Q = cbind(1, Gramm)                                # design matrix
  S = rbind(0, cbind(0, Gramm))
  M = crossprod(Q) + lambda*S
  M.inv = my.inv(M)                                  # inverse of M
  gamma.hat = crossprod(M.inv, crossprod(Q, y))
  f.hat = Q %*% gamma.hat
  A = Q %*% M.inv %*% t(Q)
  tr.A = sum(diag(A))                                # trace of hat matrix
  rss = crossprod(y - f.hat)                         # residual sum of squares
  gcv = n*rss / (n - tr.A)^2                         # obtain GCV score
  return(list(f.hat = f.hat, gamma.hat = gamma.hat, gcv = gcv))
}
lambda = 1e-8
reps = 40
lambda = 1.5^(1:reps-1) * lambda
V = sapply(lambda, function(lam) ridge.regression(X, y, lam)$gcv)
plot(lambda, V, type = 'l', main = 'GCV score', lwd = 2, ylab = 'GCV', bty='n' ) 
opt.mod = ridge.regression(X, y, lambda[which.min(V)])     # fit optimal model

### finding beta.0, beta.1 and beta.2 ##########
gamma.hat = opt.mod$gamma.hat
beta.hat.0 = opt.mod$gamma.hat[1]               # intercept
beta.hat = crossprod(gamma.hat[-1,], X)         # slope and noise term coefficients
c(beta.hat.0, beta.hat)
coef(glmnet::glmnet(X, y, alpha=0, lambda=lambda[which.min(V)]))
set.seed(3)
n = 50
x = matrix(runif(n), nrow=n, ncol=1)
x.star = as.matrix(sort(x))                      # sorted x, used by plot
y = sin(2*pi*x.star) + rnorm(n, sd=0.2)
rk.1 = function(s, t) {
  return(.5*min(s, t)^2 * max(s, t) + (1/6)*min(s, t)^3)
}
smoothing.spline = function(X, y, lambda) {
  Gramm = get.gramm(X, rkfunc=rk.1)                  # Gramm matrix (nxn)
  n = length(y)
  J = cbind(1, X)                                    # matrix with a basis for the null space of the penalty; MC note:, never name anything T (True) or t (transpose)!
  Q = cbind(J, Gramm)                                # design matrix
  m = ncol(J)                                        # dimension of the null space of the penalty
  S = matrix(0, n + m, n + m)                        # initialize S
  S[(m + 1):(n + m),(m + 1):(n + m)] = Gramm         # non-zero part of S
  M = crossprod(Q) + lambda*S
  M.inv = my.inv(M)                                  # inverse of M
  gamma.hat = crossprod(M.inv, crossprod(Q, y))
  f.hat = Q %*% gamma.hat
  A = Q %*% M.inv %*% t(Q)
  tr.A = sum(diag(A))                                # trace of hat matrix
  rss = crossprod(y - f.hat)                         # residual sum of squares
  gcv = n * rss/(n - tr.A)^2                         # obtain GCV score
  return(list(f.hat = f.hat, gamma.hat = gamma.hat,gcv = gcv))
}
lambda = 1e-8
reps = 60

lambda = 1.5^(1:reps-1) * lambda
V = sapply(lambda, function(lam) smoothing.spline(x.star, y, lam)$gcv)

# Plot of GCV
plot(1:reps, V, type = 'l', main = 'GCV score', xlab = 'i', ylab = 'GCV', bty='n') 

opt.mod.2 = smoothing.spline(x.star, y, lambda[which.min(V)])       # fit optimal model

# comparison models
mgcv::gam(y ~ s(x.star))
kernlab::gausspr(y ~ x.star, kernel='splinedot')

# Graph (Cubic Spline)
plot(x.star, y, col=scales::alpha('#ff5500', .5), pch=19, xlab = 'X', 
     ylim = c(-2.5, 1.5), xlim=c(-0.1, 1.1), ylab = 'Response', 
     main = 'Cubic Spline', bty='n', col.main='gray25')

# true
lines(x.star, sin(2*pi*x.star), lty = 1, lwd = 2, col=scales::alpha('darkred', .75))
# predictions
lines(x.star, opt.mod.2$f.hat, col = 'dodgerblue')
lines(x.star, predict(mgcv::gam(y~s(x.star))), col = 'springgreen1')
lines(x.star, kernlab::gausspr(y~x.star, kernel='splinedot')@fitted, col = 'mistyrose4')

legend(-0.1,-1.5, c('True', 'Predictions', 'gam', 'gp'), bty = 'n',
       lwd = c(2,2), col = c( scales::alpha('darkred', .75), 'dodgerblue', 'springgreen1', 'mistyrose4'),
       cex=.75)

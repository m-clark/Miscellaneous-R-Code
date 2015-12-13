## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(message=F)

## ----dummycode, echo=FALSE-----------------------------------------------
z = Z = factor(c('A','A','B','B','C','C'))
Z =  model.matrix(~ Z-1)
pander::pander(data.frame(z, Z))

## ----dataSetup, cache=TRUE-----------------------------------------------
data(sleepstudy, package='lme4')
X = model.matrix(~Days, sleepstudy)
Z = model.matrix(~factor(sleepstudy$Subject)-1)
colnames(Z) =  paste0('Subject_', unique(sleepstudy$Subject))  # for cleaner presentation later
rownames(Z) = paste0('Subject_', sleepstudy$Subject)
y = sleepstudy$Reaction

## ----mlfunc, cache=TRUE--------------------------------------------------
llMixed = function(y, X, Z, theta){
  tau = exp(theta[1])
  sigma = exp(theta[2])
  n = length(y)
  
  # evaluate covariance matrix for y
  e = tcrossprod(Z)*tau^2 + diag(n)*sigma^2
  L = chol(e)  # L'L = e
  
  # transform dependent linear model to independent
  y = backsolve(L, y, transpose=TRUE)
  X = backsolve(L, X, transpose=TRUE)
  b = coef(lm(y~X-1))
  LP = X %*% b
  
  ll = -n/2*log(2*pi) -sum(log(diag(L))) - crossprod(y-LP)/2
  -ll
}

## ----mlfuncMV, cache=TRUE------------------------------------------------
llMixedMV = function(y, X, Z, theta){
  tau = exp(theta[1])
  sigma = exp(theta[2])
  n = length(y)
  
  # evaluate covariance matrix for y
  e = tcrossprod(Z)*tau^2 + diag(n)*sigma^2

  b = coef(lm.fit(X, y))
  mu = X %*% b

  ll = -mvtnorm::dmvnorm(y, mu, e, log=T)
}

## ----optim, cache=TRUE---------------------------------------------------
paramInit = c(0, 0)
names(paramInit) = c('tau', 'sigma')

modelResults = optim(llMixed, X=X, y=y, Z=Z, par=paramInit, control=list(reltol=1e-10))
modelResultsMV = optim(llMixedMV, X=X, y=y, Z=Z, par=paramInit, control=list(reltol=1e-10))

rbind(c(exp(modelResults$par), negLogLik = modelResults$value, coef(lm(y~X-1))),
      c(exp(modelResultsMV$par), negLogLik = modelResultsMV$value, coef(lm(y~X-1)))) %>% 
  round(2)

## ----lme, cache=TRUE-----------------------------------------------------
library(lme4)
lmeMod = lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
lmeMod

## ----estRanEf, cache=TRUE------------------------------------------------
tau = exp(modelResults$par)[1]
tausq = tau^2
sigma = exp(modelResults$par)[2]
sigmasq = sigma^2
Sigma = tcrossprod(Z)*tausq/sigmasq + diag(length(y))
ranefEstimated = tausq*t(Z)%*%solve(Sigma) %*% resid(lm(y~X-1))/sigmasq
data.frame(ranefEstimated, lme4 = ranef(lmeMod)$Subject[[1]]) %>% round(2)

## ----gamData, cache=TRUE-------------------------------------------------
size = c(1.42,1.58,1.78,1.99,1.99,1.99,2.13,2.13,2.13,2.32,2.32,2.32,2.32,2.32,2.43,2.43,2.78,2.98,2.98)
wear = c(4.0,4.2,2.5,2.6,2.8,2.4,3.2,2.4,2.6,4.8,2.9,3.8,3.0,2.7,3.1,3.3,3.0,2.8,1.7)

x = size - min(size)
x = x / max(x)
d = data.frame(wear, x)

## ----gamFuncs, cache=TRUE------------------------------------------------
# function for cubic spline on [0,1]; requires data and points within domain for knots
rk <- function(x, z) {
  ((z-0.5)^2 - 1/12) * ((x-0.5)^2 - 1/12) / 4 -
    ((abs(x-z)-0.5)^4 - (abs(x-z)-0.5)^2/2 + 7/240) / 24
}

# create the model matrix
splineX <- function(x, knots) {
  q <- length(knots) + 2                    # number of parameters
  n <- length(x)                            # number of observations
  X <- matrix(1, n, q)                      # initialized model matrix
  X[, 2] <- x                               # set second column to x
  X[, 3:q] <- outer(x, knots, FUN = rk)     # remaining to cubic spline
  X
}

# set up the regression spline penalty matrix, given knot sequence knots
Sfunc = function(knots){
  q = length(knots)+2
  S = matrix(0, q, q)                       # initialize
  S[3:q, 3:q] = outer(knots, knots, FUN=rk)
  S
}

# Matrix sqrt function
matSqrt = function(S){
  UDU = eigen(S, symmetric=TRUE)
  U = UDU$vectors
  D = diag(UDU$values)
  B = crossprod(U) %*% sqrt(D)
  B
}

# the fitting function with lambda smoothing parameter
prsFit <- function(y, x, knots, lambda) {
  q = length(knots) + 2                     # dimension of basis
  n = length(x)                             # number of observations
  Xa = rbind(splineX(x, knots), 
             matSqrt(Sfunc(knots)) * sqrt(lambda))  # augmented model matrix
  y[(n + 1):(n + q)] = 0                    # augment the data vector
  lm(y ~ Xa - 1)                            # fit and return penalized regression spline
}

## ----fitCSgam, cache=TRUE, fig.align='center', out.width='50%'-----------
x_knots = 1:7/8                             # choose some knots
mod = prsFit(y=wear, x=x, knots=x_knots, lambda=.0001) # fit the penalized spline

x_pred = 0:100/100                          # values for prediction
Xp = splineX(x_pred, x_knots)               # create design matrix
predictions = Xp %*% coef(mod)

plot(x, wear, xlab='Scaled Engine size', ylab='Wear Index', pch=19,
     col="#dd4814", cex=.75, col.axis='gray50', bty='n')
lines(x_pred, predictions, col='#1e90ff')

## ----gam2mixed, cache=TRUE-----------------------------------------------
S = Sfunc(xk)
init = eigen(S)
U = init$vectors
D = diag(init$values)
poseigen = which(diag(D) > 0)  
Dpos = D[poseigen, poseigen]           # smallest submatrix containing all positive values
Xf = splineX(x, knots = xk)            # spline model matrix
U_F = U[, (ncol(U)-1):ncol(U)]         # partition eigenvector matrix
U_R = U[, 1:(ncol(U)-ncol(U_F))]
X_F = Xf %*% U_F                       # fixed part  with B_F coef to be estimated (not penalized)
X_R = Xf %*% U_R                       # random part with B_R random effects
Z = X_R %*% sqrt(Dpos)

## ----gamSleepStudy, message=FALSE, cache=TRUE----------------------------
library(mgcv); library(gamm4)
modGam = gamm4(Reaction ~ Days, random=~(1|Subject), data=sleepstudy)
summary(modGam$mer)

## ----gammSleepStudy, cache=TRUE, fig.align='center', out.width='50%', results='hold', fig.show='hold'----
modGamS = gamm4(Reaction ~ s(Days, bs='cs'), random=~(1|Subject), data=sleepstudy)
summary(modGamS$mer)    
# summary(modGamS$gam)
plot(modGamS$gam)

## ----diag, echo=FALSE, fig.align='center'--------------------------------
pal = RColorBrewer::brewer.pal(4, 'Accent')

DiagrammeR::grViz("
digraph boxes_and_circles {

  # graph statement
  graph [layout = dot,
         rankdir = TB]

  # node statements
  node [fontname = Helvetica,
        fontcolor = white,
        width = .5,
        style = filled]

  node [shape = circle,
        color = '#ff5503'] // sets as circles
  STAR 

  node [shape = box,
        color = '#1f78b4']
  Others

  node [color = '#e31a1c']
  GAMM Mixed, GAM

  node [color = '#33a02c']
  GLM, SLiM


  subgraph {
  rank = same; Mixed; GAM;
  }

  # edge statements
  STAR -> GAMM           [arrowhead=dot, color=gray]
  STAR -> Others         [arrowhead=dot, color=gray]
  GAMM -> Mixed          [arrowhead=dot, color=gray]
  GAMM -> GAM            [arrowhead=dot, color=gray]
  GAM -> GLM             [arrowhead=dot, color=gray]
  Mixed -> GLM           [arrowhead=dot, color=gray]
  GLM -> SLiM            [arrowhead=dot, color=gray]

  GAM -> Mixed           [arrowhead=none, color=gray]

}
", width=900)


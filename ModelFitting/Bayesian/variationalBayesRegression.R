
vbreg = function(X, y, a0=10e-2, b0=10e-4, c0=10e-2, d0=10e-4, tol=1e-8, maxiter=1000, ard=F){
  # X: model matrix
  # y: the response
  # a0, b0 prior parameters for tau
  # c0, d0 hyperprior parameters for alpha
  # tol: tolerance value to end iterations
  # maxiter: alternative way to end iterations
  
  
  # initializations
  X = cbind(1, X)
  D = ncol(X)
  N = nrow(X)
  w = rep(0, D) 
  XX = crossprod(X)
  Xy = crossprod(X,y)

  a_N = a0 + N/2
  
  if(!ard){
    c_N = c0 + D/2
    E_alpha = c0/d0  
  } else {
    c_N = c0 + 1/2
    E_alpha = rep(c0/d0, D)
  }
  

  tolCurrent = 1
  iter = 0
  LQ = 0
  
  while(iter < maxiter && tolCurrent > tol ){
    iter = iter + 1
    # wold = w
    
    if(!ard){
      b_N = b0 + 1/2*(crossprod(y - X%*%w) + E_alpha * crossprod(w))
      VInv = diag(E_alpha, D) + XX
      V = solve(VInv)
      w = V %*% Xy
      E_wtau = a_N/b_N*crossprod(w) + sum(diag(V))
      d_N = d0 + 1/2*E_wtau
      E_alpha = c(c_N/d_N)
    } else {
      b_N = b0 + 1/2*(crossprod(y - X%*%w) + t(w) %*% diag(E_alpha) %*% w)
      VInv = diag(E_alpha) + XX
      V = solve(VInv)
      w = V %*% Xy
      E_wtau = a_N/b_N*crossprod(w) + sum(diag(V))
      d_N = d0 + 1/2*(c(w)^2 * a_N/b_N + diag(V))
      E_alpha = c(c_N/d_N)
    }

    
    LQ_old = LQ
    suppressWarnings({
    LQ = -N/2*log(2*pi) - 1/2*(a_N/b_N * crossprod(y- crossprod(t(X), w)) + sum(XX * V)) 
    + 1/2 * log(det(V)) + D/2
    - log(gamma(a0)) + a0*log(b0) - b0*a_N/b_N + log(gamma(a_N)) - a_N*log(b_N) + a_N
    - log(gamma(c0)) + c0*log(d0) + log(gamma(c_N)) - sum(c_N*log(d_N))
    })
    tolCurrent = abs(LQ - LQ_old)
    # alternate tolerance, comment out LQ_old up to this line if using
    # tolCurrent = sum(abs(w - wold))  
  }
  
  res = list(coef=w, sigma=sqrt(1/(E_wtau/crossprod(w))), LQ=LQ, iterations=iter, tol=tolCurrent)
  if (iter>=maxiter) {
    res = append(res, warning('Maximum iterations reached.'))
  } else {res}
}


### Data set up
set.seed(1234)
n = 100
d = 3
coefs = c(1, 2, 3, 5)
sigma = 2

X = replicate(d, rnorm(n))
colnames(X) = paste0('X', 1:d)
y = cbind(1, X) %*% coefs  + rnorm(n, sd=sigma)

### Run 
res = vbreg(X, y, tol=1e-8, ard = F)
res

# With automatic relevance determination
res = vbreg(X, y, tol=1e-8, ard = T)
res


n = 150
ntest = 50
d = 100
coefs = rnorm(d + 1)
sigma = 1

Xtrain = cbind(1, replicate(d, rnorm(n)))
ytrain = Xtrain %*% coefs + rnorm(n, sd=sigma)

Xtest = cbind(1, replicate(d, rnorm(ntest)))
ytest = Xtest %*% coefs + rnorm(ntest, sd=sigma)

vbResult = vbreg(Xtrain[,-1], ytrain)
glmResult = glm.fit(Xtrain, ytrain)

vbTrainError = mean((ytrain - Xtrain %*% vbResult[['coef']])^2)
vbTestError = mean((ytest - Xtest %*% vbResult[['coef']])^2)

glmTrainError = mean((ytrain - Xtrain %*% glmResult[['coefficients']])^2)
glmTestError = mean((ytest - Xtest %*% glmResult[['coefficients']])^2)

mseResults = data.frame(vb=rbind(vbTrainError, vbTestError), 
                        glm=rbind(glmTrainError, glmTestError))
rownames(mseResults) = c('train', 'test')
mseResults


library(ggvis)

# create coefficient data set for plotting
gcoef = data.frame(wGLM=coef(glmResult), 
                   wVB=vbResult$coef, 
                   w=coefs)
gcoef = reshape2::melt(gcoef, 'w')

# same for predictions
gpred = data.frame(predGLM = Xtest %*% coef(glmResult),
                   predVB = Xtest %*% vbResult$coef,
                   y = ytest)
gpred = reshape2::melt(gpred, 'y')

gcoef %>%
  ggvis(~w, ~value) %>%
  layer_lines(~w, ~w, strokeOpacity:=.5) %>%
  layer_points(fill=~variable, fillOpacity:=.5)

gpred %>%
  ggvis(~y, ~value) %>%
  layer_lines(~y, ~y, strokeOpacity:=.5) %>%
  layer_points(fill=~variable, fillOpacity:=.75)


n = 500
ntest = 50
d = 1000
deff = 100
coefs = rnorm(deff+1)
sigma = 1

Xtrain = cbind(1, replicate(d, rnorm(n)))
ytrain = Xtrain %*% c(coefs, rep(0, d-deff))  + rnorm(n, sd=sigma)

Xtest = cbind(1, replicate(d, rnorm(ntest)))
ytest = Xtest %*% c(coefs, rep(0, d-deff))  + rnorm(ntest, sd=sigma)

vbResult = vbreg(Xtrain[,-1], ytrain, maxiter=2000)
vbARDResult = vbreg(Xtrain[,-1], ytrain, ard=T)
glmResult = glm(ytrain~., data=as.data.frame(Xtrain[,-1]))

vbTrainError = mean((ytrain - Xtrain %*% vbResult[['coef']])^2)
vbTestError = mean((ytest - Xtest %*% vbResult[['coef']])^2)

vbARDTrainError = mean((ytrain - Xtrain %*% vbARDResult[['coef']])^2)
vbARDTestError = mean((ytest - Xtest %*% vbARDResult[['coef']])^2)

glmTrainError = mean((ytrain - predict(glmResult))^2)
glmTestError = mean((ytest - predict(glmResult, as.data.frame(Xtest)))^2)

mseResults = data.frame(vb=rbind(vbTrainError, vbTestError), 
                        vbARD=rbind(vbARDTrainError, vbARDTestError), 
                        glm=rbind(glmTrainError, glmTestError))
rownames(mseResults) = c('train', 'test')
mseResults


psych::describe(vbARDResult$coef[(deff+1):d])


library(ggvis)
gcoef = data.frame(wGLM=coef(glmResult), 
                   wVB=vbResult$coef, 
                   wVBARd = vbARDResult$coef, 
                   w=c(coefs, rep(0,d-deff)))
gcoef = reshape2::melt(gcoef, 'w')

gpred = data.frame(predGLM = predict(glmResult, as.data.frame(Xtest)),
                   predVB = Xtest %*% vbResult$coef,
                   predVBARD = Xtest %*% vbARDResult$coef,
                   y = ytest)
gpred = reshape2::melt(gpred, 'y')

gcoef %>%
  ggvis(~w, ~value) %>%
  layer_lines(~w, ~w, strokeOpacity:=.5) %>%
  layer_points(fill=~variable, fillOpacity:=.5)

gpred %>%
  ggvis(~y, ~value) %>%
  layer_lines(~y, ~y, strokeOpacity:=.5) %>%
  layer_points(fill=~variable, fillOpacity:=.75)


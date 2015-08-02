# create factor structure and reproduce correlation matrix ----------------

### setup
# create initial cor/covmat
covmat = diag(1,6)
covmat[lower.tri(covmat)] = c(.6,.6,0, 0,0,
                              .6,0,0,0,
                              0,0,0,
                              .6,.6,.6)

covmat = covmat + t(covmat) - diag(diag(covmat))

# create data
d = MASS::mvrnorm(100, mu=rep(0,6), Sigma=covmat, empirical=T)

# check
round(cor(d),1)

### factor analysis
library(psych)
fact = psych::fa(d, 2)

covmatReprod = tcrossprod(loadings(fact)) + diag(fact$uniquenesses)

# compare
round(covmatReprod, 1)


### with added noise
d = MASS::mvrnorm(100, mu=rep(0,6), Sigma=covmat, empirical=F)

round(cor(d),3)

fact = psych::fa(d, 2)
covmatReprod = tcrossprod(loadings(fact)) + diag(fact$uniquenesses)

round(covmatReprod, 3)

### example of the operation to produce the cov of d[1:2,]
vecdiff = d[1,] - d[2,]
lambda = loadings(fact)
exp(-.5 * t(vecdiff) %*% (crossprod(t(lambda)) + diag(fact$uniquenesses)) %*% vecdiff)


# sim ------------------------------------------------------------

# following rasmussen code
library(pracma) # for matlab like functionality
set.seed(34)
N0 = 61^2
X = meshgrid(seq(-3,3,.1))
x = matrix(0, nrow=N0, ncol=2)
x[,1] = X[[1]]
x[,2] = X[[2]]

# Base R approach
# x = expand.grid(seq(-3,3,.1), seq(-3,3,.1))

head(x)
N = nrow(x)
L = c(1,-1) #matrix(c(1,-1,-1,1), ncol=2)  # default loadings cancel each other out
xL = x %*% L

Q = matrix(0, N, N)
Q = Q + (repmat(xL[,1,drop=F], 1, N) - repmat(t(xL[,1]), N, 1))^2
Q = Q + (repmat(x[,1,drop=F], 1, N) - repmat(t(x[,1]), N, 1))^2/36
Q = Q + (repmat(x[,2,drop=F], 1, N) - repmat(t(x[,2]), N, 1))^2/36 # 36 from l length scale set to 6^-2

# Base R approach
# Q = matrix(0, N, N)
# Q = Q + (replicate(N, xL[,1]) - t(replicate(N, xL[,1])))^2
# Q = Q + (replicate(N, x[,1]) - t(replicate(N, x[,1])))^2/36
# Q = Q + (replicate(N, x[,2]) - t(replicate(N, x[,2])))^2/36  
Q = exp(-0.5*Q)


set.seed(34)
y = t(chol(Q+1e-9*diag(N))) %*% rnorm(N)
head(y)
# y = MASS::mvrnorm(N, mu=rep(0,N), Sigma= Q+1e-9*diag(N))


plot(x[,1],y)
plot(x[,2],y)
library(ggplot2); library(scales)
pal = scale_color_gradient2(low = muted("red"), mid = "white",high = muted("blue"), midpoint=0)
car::scatter3d(y ~ x[,1] + x[,2], point.col=rainbow(N), surface=F)
car::scatter3d(y ~ x[,1] + x[,2], point.col=pal$palette((y+abs(min(y)))/max(y+abs(min(y)))), surface=F)


# original matlab code
# n = 61^2; D=2;
# [x1,x2]=meshgrid(-3:0.1:3,-3:0.1:3);
# x=[x1(:),x2(:)];
# L = [1 -1];
# xL = x*L';
# Q = zeros(n,n);
# Q = Q + (repmat(xL(:,1),1,n)-repmat(xL(:,1)',n,1)).^2;
# Q = Q + (repmat(x(:,1),1,n)-repmat(x(:,1)',n,1)).^2/36;
# Q = Q + (repmat(x(:,2),1,n)-repmat(x(:,2)',n,1)).^2/36;
# Q = exp(-0.5*Q);
# randn('seed',34)
# y = chol(Q+1e-9*eye(n))'*randn(n,1);


# run stan ----------------------------------------------------------------

# redo for more manageable data; In addition, the above factor analytic model is
# unidentified unless you fix one of L in the estimation process (and
# technically you need another x also to get beyond 'just-identified'); the stan
# model code does fix the first value to 1, in keeping with standard structural
# equation modeling approaches. If you don't, you'll see that at least one or
# two chains start to flip flop assigning the loading estimates, but otherwise
# do ok.

set.seed(34)

N0 = 15
x = as.matrix(expand.grid(seq(-3,3, length=N0), seq(-3,3, length=N0)))
head(x)
N = nrow(x)
L = matrix(c(1,-1), ncol=1)      # this model is unidentified from an fa perspective
xL = x %*% L

Q = matrix(0, N, N)
Q = Q + (replicate(N, xL[,1]) - t(replicate(N, xL[,1])))^2
Q = Q + (replicate(N, x[,1]) - t(replicate(N, x[,1])))^2/36
Q = Q + (replicate(N, x[,2]) - t(replicate(N, x[,2])))^2/36  # 36 from l length scale set to 6
Q = exp(-0.5*Q)

sigma = .01  # adding a little more noise so stan doesn't have to estimate a boundary (it actually did though, just made it really slow)
y = t(chol(Q+sigma*diag(N))) %*% rnorm(N)
head(y)
psych::describe(y)

plot(x[,1],y)
plot(x[,2],y)
car::scatter3d(x[,1], y, x[,2], surface=F, point.col=rainbow(N))

# xtest = seq(x)

standata = list(N=nrow(x), D=ncol(x), X=x, y=y[,1], K=1, Xtest=x, Ntest=nrow(x))

library(rstan)
fit0 = stan(data=standata, file = 'ModelFitting/gp Examples/gpStan_squaredExponentialFactorAnalysis.stan', iter = 1, chains=1)

# timed fit for 1000 iterations(12 min for N0=15, 2.5 for N0=10)
# p = proc.time()
# fit = stan(data=standata, file = 'ModelFitting/gp Examples/gpStan_squaredExponentialFactorAnalysis.stan', iter = 1000, warmup = 200, chains=1, fit = fit0)
# (proc.time() - p)/60

# parallelize chains
iterations = 6000
wu = 1000
th = 20
chains = 4


p = proc.time()
fit = stan(data=standata, iter = iterations, warmup = wu, thin=th, chains=chains, fit = fit0, cores=chains)
(proc.time() - p)/3600



# Summarize and Vis -------------------------------------------------------

# takes a bit to print
print(fit, digits=3, par=c('eta_sq','sigma_sq','l_sq','lambda'))
traceplot(fit, par=c('eta_sq','sigma_sq','l_sq','lambda'), inc_=F)

# library(shinyStan)
# launch_shinystan(fit)

# Visualize fit
yRepMean = get_posterior_mean(fit, 'yRep')[,5]
gdat = data.frame(x, y=yRepMean)
car::scatter3d(y ~., data=gdat, point.col=rainbow(N), surface=F)



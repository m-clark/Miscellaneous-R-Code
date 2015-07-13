# create factor structure and reproduce correlation matrix ----------------

### setup
covmat = diag(1,6)
covmat[lower.tri(covmat)] = c(.6,.6,0, 0,0,
                              .6,0,0,0,
                              0,0,0,
                              .6,.6,.6)
covmat[upper.tri(covmat)] = t(covmat)[upper.tri(t(covmat))]

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

# example of the operation to produce the cov of d[1:2,]
vecdiff = d[1,] - d[2,]
lambda = loadings(fact)
exp(-.5 * t(vecdiff) %*% (crossprod(t(lambda)) + diag(fact$uniquenesses)) %*% vecdiff)


# sim ------------------------------------------------------------

# following rasmussen code
library(pracma) # for matlab like functionality
set.seed(34)
N =
X = meshgrid(seq(-3,3,.1))
x = matrix(0, nrow=length(seq(-3,3,.1))^2, ncol=2)
x[,1] = X[[1]]
x[,2] = X[[2]]
head(x)
N = nrow(x)
L = c(1,-1) #matrix(c(1,-1,-1,1), ncol=2)  # default loadings cancel each other out
xL = x %*% L

Q = matrix(0, N, N)
Q = Q + (replicate(N, xL[,1]) - t(replicate(N, xL[,1])))^2
Q = Q + (replicate(N, x[,1]) - t(replicate(N, x[,1])))^2/36
Q = Q + (replicate(N, x[,2]) - t(replicate(N, x[,2])))^2/36  # 36 from l length scale set to 6^-2
Q = exp(-0.5*Q)


set.seed(34)
y = t(chol(Q+1e-9*diag(N))) %*% rnorm(N)
head(y)
# y = MASS::mvrnorm(N, mu=rep(0,N), Sigma= Q+1e-9*diag(N))


plot(x[,1],y)
plot(x[,2],y)
car::scatter3d(x[,1], y, x[,2], surface=F, point.col=rainbow(N))

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
X = expand.grid(seq(-3,3, length=N0), seq(-3,3, length=N0))  # base r
x = matrix(0, nrow=N0^2, ncol=2)
x[,1] = X[[1]]
x[,2] = X[[2]]
head(x)
N = nrow(x)
L = c(1,-1)      # this model is unidentified from an fa perspective
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

xtest = x

standata = list(N=nrow(x), D=ncol(x), X=x, y=y[,1], K=1, Xtest=xtest, Ntest=nrow(xtest))

library(rstan)
fit0 = stan(data=standata, file = 'gpStan_squaredExponentialFactorAnalysis.stan', iter = 1, chains=1)


# parallelize chains
iterations = 12000
wu = 2000
th = 40
chains = 4

library(parallel)
cl = makeCluster(4)
clusterExport(cl, c('standata', 'fit0','iterations', 'wu', 'th', 'chains'))
clusterEvalQ(cl, library(rstan))


p = proc.time()
fit = parLapply(cl, seq(chains), function(chain) stan(data=standata,
                                                      iter = iterations,
                                                      warmup = wu, thin=th,
                                                      chains=1, chain_id=chain,
                                                      fit = fit0)
)
(proc.time() - p)/3600

stopCluster(cl)



# Summarize and Vis -------------------------------------------------------

fit = sflist2stanfit(fit)
# takes a bit to print
print(fit, digits=3)
# traceplot(fit)

library(shinyStan)
launch_shinystan(fit)


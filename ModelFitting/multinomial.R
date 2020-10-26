# for more detail on this see my categorical document at m-clark.github.io/docs/logregmodels.html

# Import data and setup  --------------------------------------------------


library(haven)
library(tidyverse)
library(mlogit)

program = read_dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta") %>%
  as_factor() %>%
  mutate(prog = relevel(prog, ref = "academic"))

head(program[, 1:5])


# convert to long form for mlogit
programLong = program %>%
  select(id, prog, ses, write) %>%
  mlogit.data(
    data = ,
    shape = 'wide',
    choice = 'prog',
    id.var = 'id'
  )

head(programLong)

mlogit_mod = mlogit(prog ~ 1| write + ses, data = programLong)
mlogit_coefs = coef(mlogit_mod)[c(1,5,7,3,2,6,8,4)]



# A basic loglik approach for comparison ----------------------------------

multinomregML <- function(par, X, y) {
  levs = levels(y)
  ref  = levs[1]              # reference level (category label 1)
  
  y0 = y == ref
  y1 = y == levs[2]             # category 2
  y2 = y == levs[3]             # category 3
  
  beta = matrix(par, ncol = 2)
  
  # more like mnlogit package depiction in its function
  # V1 = X %*% beta[ ,1]
  # V2 = X %*% beta[ ,2]
  # ll = sum(-log(1 + exp(V1) + exp(V2))) + sum(V1[y1], V2[y2])
  
  V = X %*% beta                           # a vectorized approach
  baseProbVec = 1 / (1 + rowSums(exp(V)))  # reference group probabilities
  
  loglik = sum(log(baseProbVec))  + crossprod(c(V), c(y1, y2))
  
  loglik
}


out = optim(
  runif(8,-.1, .1),
  multinomregML,
  X = model.matrix(prog ~ ses + write, data = program),
  y = program$prog,
  control = list(
    maxit   = 1000,
    reltol  = 1e-12,
    ndeps   = rep(1e-8, 8),
    trace   = TRUE,
    fnscale = -1,
    type    = 3
  ),
  method = 'BFGS'
)

# out$par

cbind(out$par, mlogit_coefs) %>% 
  round(4)

# setup for loglike comparison
X = model.matrix(prog ~ ses + write, data = program)
y = program$prog
pars = matrix(out$par, ncol = 2)
V = X %*% pars
acadprob = 1 / (1+rowSums(exp(V)))
fitnonacad = exp(V) * matrix(rep(acadprob, 2), ncol = 2)
fits = cbind(acadprob, fitnonacad)
yind = model.matrix( ~ -1 + prog, data = program)


# because dmultinom can't take matrix for prob
ll = 0

for (i in 1:200){
  ll = ll + dmultinom(yind[i, ],
                      size = 1,
                      prob = fits[i, ],
                      log  = TRUE)
}

ll

out$value

logLik(mlogit_mod)



# Alternative specific and constant variables -----------------------------

# in this example, price is alternative invariant (Z) income is
# individual/alternative specific (X), and catch is alternative specific (Y)

library(mnlogit)
data(Fish)

head(Fish)

fm  = formula(mode ~ price | income | catch)

fit = mnlogit(fm, Fish)
# fit = mlogit(fm, Fish)
summary(fit)


# X dim nrow(Fish)/K x p + 1 (intercept)
# Z, Y nrow(N); Y has alt specific coefs; then for Z ref group dropped so nrow = nrow*(K-1)/K
# for ll everything through previous X the same
# then calc probmat for Y and Z, add to X probmat, and add to base

multinomregML2 <- function(par, X, Y, Z, respVec, choice) {
  
  N = sum(choice)
  K = length(unique(respVec))
  levs = levels(respVec)
  
  xpar = matrix(par[1:6], ncol = K-1)
  ypar = matrix(par[7:10], ncol = K)
  zpar = matrix(par[length(par)], ncol = 1)
  
  # Calc X
  Vx  = X %*% xpar
  
  # Calc Y (mnlogit finds N x 1 results by going through 1:N, N+1:N*2 etc; then
  # makes 1 vector, then subtracts the first 1:N from whole vector, then makes
  # Nxk-1 matrix with N+1:end values (as 1:N are just zero)); creating the
  # vector and rebuilding the matrix is unnecessary though
  Vy = sapply(1:K, function(alt) 
    Y[respVec == levs[alt], , drop = FALSE] %*% ypar[alt])
  
  Vy = Vy[,-1] - Vy[,1]
  
  # Calc Z
  Vz = Z %*% zpar
  Vz = matrix(Vz, ncol = 3)
  
  # all Vs must fit into N x K -1 matrix where N is nobs (i.e. individuals)
  V = Vx + Vy + Vz
  
  ll0 = crossprod(c(V), choice[-(1:N)])
  baseProbVec <- 1 / (1 + rowSums(exp(V)))
  loglik = sum(log(baseProbVec)) + ll0
  
  loglik
  
  # note fitted values via
  # fitnonref = exp(V) * matrix(rep(baseProbVec, K-1), ncol = K-1)
  # fitref = 1-rowSums(fitnonref)
  # fits = cbind(fitref, fitnonref)
}


inits = runif(11, -.1, .1)
mdat  = mnlogit(fm, Fish)$model  # this data already ordered!

# X has a constant value across alternatives; the coefficients regard selection of alternative relative to reference
X = cbind(1, mdat[mdat$`_Alt_Indx_` == 'beach', 'income'])
dim(X)
head(X)

# Y will use the complete data to start; coefficients will be differences from reference alternative coefficient
Y = as.matrix(mdat[, 'catch', drop = FALSE])
dim(Y)

# Z are difference scores from reference group
Z = as.matrix(mdat[mdat$`_Alt_Indx_` != 'beach', 'price', drop = FALSE])
Z = Z - mdat[mdat$`_Alt_Indx_` == 'beach', 'price']
dim(Z)

respVec = mdat$`_Alt_Indx_` # first 10 should be 0 0 1 0 1 0 0 0 1 1 after beach dropped

# debugonce(multinomregML2)
multinomregML2(inits, X, Y, Z, respVec, choice = mdat$mode)

out = optim(
  par = rep(0, 11),
  multinomregML2,
  X = X,
  Y = Y,
  Z = Z,
  respVec = respVec,
  choice  = mdat$mode,
  control = list(
    maxit   = 1000,
    reltol  = 1e-12,
    ndeps   = rep(1e-8, 11),
    trace   = TRUE,
    fnscale = -1,
    type    = 3
  ),
  method = 'BFGS'
)

# out

# round(out$par, 3)

round(cbind(out$par, coef(fit)), 3)

cbind(logLik(fit), out$value)

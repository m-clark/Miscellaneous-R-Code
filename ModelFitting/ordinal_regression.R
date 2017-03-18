# The following demonstrates a standard cumulative link ordinal regression model
# via maximum likelihood. Default is with probit link function. Alternatively
# you can compare it with a logit link, which will result in values roughly 
# 1.7*parameters estimates from the probit.


ll_ord = function(par, X, y, probit=T) {
  K = length(unique(y))                        # number of classes K
  ncuts = K-1                                  # number of cutpoints/thresholds
  cuts = par[(1:ncuts)]                        # cutpoints
  beta = par[-(1:ncuts)]                       # regression coefficients
  lp = X %*% beta                              # linear predictor
  ll = rep(0, length(y))                       # log likelihood
  pfun = ifelse(probit, pnorm, plogis)         # which link to use
  
  for(k in 1:K){
    if (k==1) {
      ll[y==k] = pfun((cuts[k] - lp[y==k]), log = TRUE)
    }
    else if (k < K) {
      ll[y==k] = log(pfun(cuts[k] - lp[y==k]) - pfun(cuts[k-1] - lp[y==k]))
    }
    else {
      ll[y==k] = log(1 - pfun(cuts[k-1] - lp[y==k])) 
    }
  }
  return(-sum(ll))
}

# data generation from the probit perspective, where the underlying continuous
# latent variable is normally distributed

set.seed(808)
N = 1000                                       # Sample size
x = cbind(x1 = rnorm(N), x2 = rnorm(N))        # predictor variables
beta = c(1,-1)                                 # coefficients
y_star = rnorm(N, mean=x %*% beta)             # the underlying latent variable
y_1 = y_star > -1.5                            # -1.5 first cutpoint
y_2 = y_star > .75                             # .75 second cutpoint
y_3 = y_star > 1.75                            # 1.75 third cutpoint
y = y_1 + y_2 + y_3 + 1                        # target

table(y)

d = data.frame(x, y=factor(y))

init = c(-1,1,2, 0,0)                          # initial values

result_probit = optim(init, ll_ord, y = y, X = x, probit=T,
                      control=list(reltol=1e-10))
result_logit  = optim(init, ll_ord, y = y, X = x, probit=F,
                      control=list(reltol=1e-10))

# compare with ordinal package
library(ordinal)
result_ordpack_probit = clm(y ~ x1 + x2, data=d, link='probit')
result_ordpack_logit  = clm(y ~ x1 + x2, data=d, link='logit')


resprobit = data.frame(method=c('ll_ord', 'ordpack'),
                       rbind(coef(result_probit), coef(result_ordpack_probit)))
colnames(resprobit) = c('method', paste0('cut_', 1:3), 'beta1', 'beta2')
pander::pander(resprobit, round=3)

reslogit = data.frame(method=c('ll_ord', 'ordpack'),
                       rbind(coef(result_logit), coef(result_ordpack_logit)))
colnames(reslogit) = c('method', paste0('cut_', 1:3), 'beta1', 'beta2')
pander::pander(reslogit, round=3)

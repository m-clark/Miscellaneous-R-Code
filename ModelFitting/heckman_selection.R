# heckman selection model based on bleven's example here: 
# https://www3.nd.edu/~wevans1/ecoe60303/sample_selection_example.ppt, which is 
# more or less the 'classic' example, variations of which you'll find all over
# regarding women's wages.
#
# Draw 10,000 obs at random 
# educ uniform over [0,16] 
# age uniform over [18,64]
# wearnl=4.49 + 0.08*educ + 0.012*age + ε 
# 
# Generate missing data for wearnl drawn
# z from standard normal [0,1]  z is actually never explained in the slides, I
# think it's left out on slide 3 and just represents an additional covariate
# 
# d*=-1.5+0.15*educ+0.01*age+0.15*z+v wearnl missing if d*≤0 wearn reported if 
# d*>0
# 
# wearnl_all=wearnl with non-missing obs



# Data Setup --------------------------------------------------------------

set.seed(123456)
N = 10000
educ = sample(1:16, N, replace = T)
age = sample(18:64, N, replace = T)

covmat = matrix(c(.46^2,.25*.46, .25*.46, 1), ncol=2)
errors = mvtnorm::rmvnorm(N, sigma=covmat)
z = rnorm(N)
e = errors[,1]
v = errors[,2]
wearnl = 4.49 + .08*educ + .012*age + e
d_star = -1.5 + 0.15*educ + 0.01*age + 0.15*z + v
observed_index  = d_star > 0

d = data.frame(wearnl, educ, age, z, observed_index)



# Comparison models -------------------------------------------------------

# lm based on full data
lm_all = lm(wearnl ~ educ + age, data=d)

# lm based on observed data
lm_obs = lm(wearnl ~ educ + age, data=d[observed_index,])

summary(lm_all)
summary(lm_obs) # smaller coefs, resid standard error



# 2 step approach ---------------------------------------------------------

# The two-step approach first conducts a probit model regarding whether the
# individual is observed or not, in order to calculate the inverse mills ratio,
# or 'nonselection hazard'. The second step is a standard lm.

## Step 1: probit model
probit = glm(observed_index ~ educ + age + z, data=d, family=binomial(link = 'probit'))
summary(probit)


# http://www.stata.com/support/faqs/statistics/inverse-mills-ratio/
probit_lp = predict(probit)
mills0 = dnorm(probit_lp)/pnorm(probit_lp)
summary(mills0)
# identical formulation 
# probit_lp = -predict(probit)
# imr = dnorm(probit_lp)/(1-pnorm(probit_lp))
imr = mills0[observed_index]
summary(imr)
ggplot2::qplot(imr, geom='histogram')

## Step 2: standard regression model using the inverse mills ratio as covariate
lm_select = lm(wearnl ~ educ + age + imr, data=d[observed_index,])
summary(lm_select)


## compare to sampleSelection package
library(sampleSelection)
selection_2step = selection(observed_index ~ educ + age + z, wearnl ~ educ + age, method = '2step')
summary(selection_2step)
coef(lm_select)['imr']/summary(lm_select)$sigma                   # slightly off
coef(lm_select)['imr']/summary(selection_2step)$estimate['sigma','Estimate']



# Maximum Likelihood ------------------------------------------------------

# likelihood function takes arguments as follows:
# par: the regression coefficients pertaining to the two models, the residual se
# sigma and rho for the correlation estimate
# X: observed data model matrix for the linear regression model
# Z: complete data model matrix for the probit model
# y: the target variable
# observed_index: an index denoting whether y is observed

ll_select <- function(par, X, Z, y, observed_index) {
  gamma = par[1:4]
  lp_probit = Z %*% gamma
  
  beta = par[5:7]
  lp_lm = X %*% beta
  
  sigma = par[8]
  rho = par[9]
  
  ll = sum(log(1-pnorm(lp_probit[!observed_index]))) +
    - log(sigma) +
    sum(dnorm(y, mean = lp_lm, sd = sigma, log = T)) +
    sum(pnorm((lp_probit[observed_index] + rho/sigma * (y-lp_lm)) / sqrt(1-rho^2), log.p = T))
  - ll
}

X = model.matrix(lm_select)
Z = model.matrix(probit)

init = c(coef(probit), coef(lm_select)[-4],  1, 0)  # initial values

# without bounds for sigma and rho you'll get warnings, but does fine anyway
heckman_ml_unbounded = optim(init, ll_select, 
                             X=X[,-4], 
                             Z=Z, 
                             y=wearnl[observed_index], 
                             observed_index = observed_index,
                             method='BFGS',
                             control=list(maxit=1000, reltol=1e-15), 
                             hessian=T)
heckman_ml_bounded = optim(init, ll_select, 
                           X=X[,-4], 
                           Z=Z, 
                           y=wearnl[observed_index], 
                           observed_index = observed_index,
                           method='L-BFGS', 
                           lower=c(rep(-Inf, 7), 1e-10, -1), 
                           upper=c(rep(Inf, 8), 1),
                           control=list(maxit=1000, factr=1e-15), 
                           hessian=T)
# heckman_ml_unbounded
# heckman_ml_bounded

# comparison model
selection_ml = selection(observed_index ~ educ + age + z, wearnl ~ educ + age, method = 'ml')
# summary(selection_ml)




# Compare results ---------------------------------------------------------

library(tidyverse)

# compare coefficients
data_frame(model = rep(c('probit', 'lm', 'both'), c(4,4,1)),
           par=names(coef(selection_ml)),
           sampselpack_ml=coef(selection_ml),
           unbounded_ml=heckman_ml_unbounded$par,
           bounded_ml=heckman_ml_bounded$par,
           explicit_twostep = c(coef(probit), coef(lm_select)[1:3], summary(lm_select)$sigma, 
                                coef(lm_select)['imr']/summary(lm_select)$sigma),
           sampselpack_2step = coef(selection_2step)[-8]
           ) %>%
  mutate_if(is.numeric, round, digits=3)

# compare standard errors
data_frame(model = rep(c('probit', 'lm', 'both'), c(4,4,1)),
           par=names(coef(selection_ml)),
           sampselpack_ml=sqrt(diag(solve(-selection_ml$hessian))),
           unbounded_ml=sqrt(diag(solve(heckman_ml_unbounded$hessian))),
           bounded_ml= sqrt(diag(solve(heckman_ml_bounded$hessian))),
           explicit_twostep = c(summary(probit)$coefficients[,2], summary(lm_select)$coefficients[-4,2], NA, NA),
           sampselpack_2step = summary(selection_2step)$estimate[-8,2]
           ) %>%
  mutate_if(is.numeric, round, digits=3)

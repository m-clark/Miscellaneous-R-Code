# Demonstration of a simple marginal structural model for estimation of
# so-called 'causal' effects using inverse probability weighting.

# Example data is from, and comparison made to, the ipw package.  See more here:
# https://www.jstatsoft.org/article/view/v043i13/v43i13.pdf


# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(ipw)

# Data Setup --------------------------------------------------------------

# this example is from the helpfile at ?ipwpoint
set.seed(16)
n <- 1000
simdat <- data.frame(l = rnorm(n, 10, 5))
a.lin <- simdat$l - 10
pa <- plogis(a.lin)

simdat <- simdat %>% 
  mutate(
    a = rbinom(n, 1, prob = pa),
    y = 10 * a + 0.5 * l + rnorm(n, -10, 5)
  )


ipw_result <- ipwpoint(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ l,
  data = simdat
)

summary(ipw_result$ipw.weights)
ipwplot(ipw_result$ipw.weights)


# Create the weights by hand for demonstration ----------------------------

ps_num = fitted(glm(a ~ 1, data = simdat, family = 'binomial'))
ps_num[simdat$a == 0] = 1 - ps_num[simdat$a == 0]

ps_den = fitted(glm(a ~ l, data = simdat, family = 'binomial'))
ps_den[simdat$a == 0] = 1 - ps_den[simdat$a == 0]

wts = ps_num / ps_den

# compare
rbind(summary(wts), summary(ipw_result$ipw.weights))

# Add inverse probability weights to the data if desired
simdat <- simdat %>% 
  mutate(sw = ipw_result$ipw.weights)


# Marginal Structural Model -----------------------------------------------

# Marginal structural model for the causal effect of `a` on `y` corrected for
# confounding by `l` using inverse probability weighting with robust standard
# error from the survey package.

library("survey")

msm <- svyglm(
  y ~ a,
  design = svydesign(~ 1, weights = ~ sw, data = simdat)
)

summary(msm)

# create the likelihood function for using the weights
maxlike = function(
  par,             # parameters to be estimated; first is taken to be sigma
  X,               # model matrix
  y,               # target variable
  wts              # estimated weights
) {
  beta = par[-1]
  lp = X %*% beta
  sigma = exp(par[1])  # eponentiated value to stay positive
  ll = dnorm(y, mean = lp, sd = sigma, log = TRUE)  # weighted likelihood
  
  -sum(ll*wts)
  
  # same as
  # ll = dnorm(y, mean = lp, sd = sigma)^wts
  # -sum(log(ll))
}

X = cbind(1, simdat$a)
y = simdat$y

result = optim(
  par = c(sigma = 0, intercept = 0, b = 0),
  fn  = maxlike,
  X   = X,
  y   = y,
  wts = wts,
  hessian = TRUE,
  method  = 'BFGS',
  control = list(abstol = 1e-12)
)

dispersion = exp(result$par[1])^2
beta = result$par[-1]


# Compute standard errors -------------------------------------------------

# the following is the survey package raw version to get the appropriate
# standard errors, which the ipw approach uses
glm_basic = glm(y ~ a, data = simdat, weights = wts)     # to get unscaled cov
res = resid(glm_basic, type = 'working')                 # residuals
glm_vcov_unsc = summary(glm_basic)$cov.unscaled          # weighted vcov unscaled by dispersion solve(crossprod(qr(X)))
estfun = X * res * wts                  
x = estfun %*% glm_vcov_unsc 

# get standard errors
se = sqrt(diag(crossprod(x)*n/(n-1)))                    # a robust standard error
se_robust = sqrt(diag(sandwich::sandwich(glm_basic)))    # much easier way to get it
se_msm    = sqrt(diag(vcov(msm)))

tibble(se, se_robust, se_msm)


# Compare results ---------------------------------------------------------

tibble(
  Estimate  = beta,
  init_se   = sqrt(diag(solve(result$hessian)))[c('intercept', 'b')],   # same as scaled se from glm_basic
  se_robust = se_robust,
  t = Estimate/se,
  p = 2*pt(abs(t), df = n-ncol(X), lower.tail = FALSE),  
  dispersion = dispersion             
)

# compare to msm
broom::tidy(msm)

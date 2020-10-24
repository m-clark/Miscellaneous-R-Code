# A simple demonstration of tobit regression via maximum likelihood. The issue
# is one where data is censored such that while we observe the value, it is not
# the true value, which would extend beyond the range of the observed data. This
# is very commonly seen in cases where the dependent variable has been given
# some arbitrary cutoff at the lower or upper end of the range, often resulting
# in floor or ceiling effects respectively. The conceptual idea is that we are
# interested in modeling the underlying latent variable that would not have such
# restriction if it was actually observed.


tobit <- function(par, X, y, ul = -Inf, ll = Inf) {
  
  # this function only takes a lower OR upper limit
  
  # parameters
  sigma = exp(par[length(par)]) 
  beta  = par[-length(par)]
  
  # create indicator depending on chosen limit
  if (!is.infinite(ll)) {
    limit = ll
    indicator = y > ll
  } else {
    limit = ul
    indicator = y < ul
  }
  
  # linear predictor
  lp = X %*% beta
  
  # log likelihood
  ll = sum(indicator * log((1/sigma)*dnorm((y-lp)/sigma)) ) + 
    sum((1-indicator) * log(pnorm((lp-limit)/sigma, lower=is.infinite(ll))))
  
  -ll
}



# demonstrate censoring with an upper limit -------------------------------

# Data setup
# Data regards academic aptitude (GRE scores) with will be modeled using reading
# and math test scores, as well as the type of program the student is enrolled
# in (academic, general, or vocational).  See this for an applied example and
# more detail- https://stats.idre.ucla.edu/r/dae/tobit-models/

library(tidyverse)

acad_apt = read_csv("https://stats.idre.ucla.edu/stat/data/tobit.csv") %>%
  mutate(prog = factor(prog, labels = c('acad', 'general', 'vocational')))

# setup data and initial values

initmod = lm(apt ~ read + math + prog, data=acad_apt)
X = model.matrix(initmod)
init = c(coef(initmod), log_sigma=log(summary(initmod)$sigma))

res = optim(
  par = init,
  tobit,
  y  = acad_apt$apt,
  X  = X,
  ul = 800,
  method  = 'BFGS',
  control = list(maxit = 2000, reltol = 1e-15)
)


# this would be more akin to the default Stata default approach
# optim(par=init, tobit, y=acad_apt$apt, X=X, ul=800,
#       control=list(maxit=16000, reltol=1e-15))


# compare to AER package tobit function

library(survival) # from base R

aer_mod = AER::tobit(
  apt ~ read + math + prog,
  data = acad_apt,
  left = -Inf,
  right = 800
)

rbind(
  tobit = c(
    res$par[1:5],
    sigma = exp(res$par[6]),
    logLike = -res$value
  ),
  AER = c(coef(aer_mod), aer_mod$scale, logLik(aer_mod))
) %>% 
  round(3)

# AER is actually just using survreg from the survival package. Survival models
# are usually for modeling time to some event, e.g. death in medical studies,
# and the censoring comes from the fact that the observed event does not occur
# for some people. Like our tobit function, an indicator is needed to denote who
# is or isn't censored. In survival models, the indicator is for the event
# itself, and means they are NOT censored.  So we'll reverse the indicator used
# in the tobit function for survreg.

surv_mod = survreg(Surv(apt, apt < 800, type = 'right') ~ read + math + prog,
                   data = acad_apt,
                   dist = 'gaussian')

# Compare all results

rbind(
  tobit = c(
    res$par[1:5],
    sigma = exp(res$par[6]),
    logLike = -res$value
  ),
  AER = c(coef(aer_mod), aer_mod$scale, logLik(aer_mod)),
  survival = c(coef(surv_mod), surv_mod$scale, logLik(surv_mod))
) %>% 
  round(3)


# Demonstrate censoring with an lower limit -------------------------------

# create a censored data situation for the low end.  The scale itself would be
# censored for anyone scoring a 200, but that basically doesn't happen. In this
# data, 15 are less than a score of 500, so we'll do that.

acad_apt = acad_apt %>%
  mutate(apt2 = apt,
         apt2 = if_else(apt2 < 500, 500, apt2))

res = optim(
  par = init,
  tobit,
  y  = acad_apt$apt2,
  X  = X,
  ll = 400,
  method  = 'BFGS',
  control = list(maxit = 2000, reltol = 1e-15)
)

aer_mod = AER::tobit(apt2 ~ read + math + prog,
                     data = acad_apt,
                     left = 400)


rbind(
  tobit = c(
    res$par[1:5],
    sigma = exp(res$par[6]),
    logLike = -res$value
  ),
  AER = c(coef(aer_mod), aer_mod$scale, logLik(aer_mod))
) %>% 
  round(3)


library(tidyverse); library(rstan); library(shinystan); library(ltm)

# Data setup
data(Abortion, package='ltm')
Abortion_ = sapply(Abortion, as.numeric)
datalist = list(Y=Abortion_, N=nrow(Abortion_), J=ncol(Abortion_), Model=1)

# 1PM
test_rasch = stan(file='IRT_rasch.stan', data=datalist, iter=10)
real_rasch = stan(fit=test_rasch, data=datalist, cores=4, iter=3000, warmup=1000, thin=8)

print(real_rasch, par=c('difficulty', 'discrim'), digits=3)

# rasch(Abortion, constraint=cbind(ncol(Abortion) + 1, 1))
rasch(Abortion)

# std of random effect === discrimination; coefs are comparable to rasch(..., IRT.param=FALSE)
# Abortion %>% 
#   mutate(Subject=1:nrow(Abortion)) %>% 
#   gather(key=Item, value=Response, -Subject) %>% 
#   lme4::glmer(Response ~ -1 + Item + (1|Subject), data=., family=binomial)
# or
# brms::brm(Response ~ -1 + Item + (1|Subject), data=., family=bernoulli)


# launch_shinystan(real_rasch)


# 2PM
test_2pm = stan(file='IRT_2PM.stan', data=datalist, iter=10)
real_2pm = stan(fit=test_2pm, data=datalist, cores=4, iter=3000, warmup=1000, thin=8)

print(real_2pm, par=c('difficulty', 'discrim'), digits=3)

ltm(Abortion ~ z1)
# launch_shinystan(real_2pm)


# 3PM
test_3pm = stan(file='IRT_3PM.stan', data=datalist, iter=10)
real_3pm = stan(fit=test_3pm, data=datalist, cores=4, iter=4000, warmup=2000, thin=8)

print(real_3pm, par=c('difficulty', 'discrim', 'guess'))

tpm(Abortion)
# launch_shinystan(real_3pm)


# 4PM
test_4pm = stan(file='IRT_4PM.stan', data=datalist, iter=10)
real_4pm = stan(fit=test_4pm, data=datalist, cores=4, iter=4000, warmup=2000, thin=8)

library(sirt)
compare4pm = rasch.mml2(Abortion, est.a = 1:4 , est.c=1:4 , est.d = 1:4,
                       min.a=0, min.b=-3, max.b=3, min.d = .95 , max.c = .25)  # so nice to have no verbose option!
compare4pm$item[, c('a','b','c','d')] # a = discrim, b = diff, c = guess, d = ceiling
print(real_4pm, par=c('difficulty', 'discrim', 'guess', 'ceiling'))

# launch_shinystan(real_4pm)

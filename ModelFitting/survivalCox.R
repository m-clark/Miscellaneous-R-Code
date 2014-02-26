# Some simple demonstrations of a standard Cox, Cox with time-varying covariates 
# and a stratified Cox

####################
### Standard Cox ###
####################

### Create the Data ### 
set.seed(12)

dur = 1:10
kittyblarg= rnorm(10)                                 # something happened to kitty!
kittyhappy = rep(0:1,times=5)                         # is kitty happy?
kittydied = sample(0:1,10,replace=T)                  # kitty died! oh noes!
d = data.frame(kittyblarg,kittyhappy,dur,kittydied)

# Inspect
d


### Create a function to feed to optim ###
pl = function(pars, preds, died, t) {
  # Arguments- pars: coefficients of interest; preds: predictor matrix; 
  # died: death; t: time
  b = pars
  X = as.matrix(preds[order(t),])
  died2 = died[order(t)]
  LP = X%*%b                                                 # Linear predictor
  ll = numeric(nrow(X))                                      # initiallize log likelihood due to looping, not necessary
  rows = 1:nrow(preds)
  for (i in rows){
    riskset = ifelse(rows < i, F, T)                         # identify risk set
    ll[i] = died2[i]*(LP[i] - log(sum(exp(LP[riskset]))) )   # log likelihood
  }
  -sum(ll)
}


### Run the model via optim and compare to survival ###
init.val = c(0,0)
out = optim(par = init.val, fn = pl, preds = d[,c('kittyblarg','kittyhappy')],
            died = d[,'kittydied'], t=dur, 
            method="BFGS", hessian=T)
out
# create a summary table
B = out$par
se = sqrt(diag(solve(out$hessian)))
Z = B/se
data.frame(B, exp=exp(B), se, Z, 
           p=ifelse( Z>0, pnorm(Z, lower=F)*2, pnorm(Z, lower=T)*2) )

# Compare to survival package
library(survival)
coxmod = coxph(Surv(dur, kittydied) ~ kittyblarg + kittyhappy)
coxmod; coxmod$loglik[2]
# summary(coxmod)



#################################
### Time-varying coefficients ###
#################################
# Note that technically nothing new is going on here relative to the previous model.  
# See the vignette for the survival package for further details.

### Create Data ###
set.seed(123)

# In the following we'll first create some noisy time points
t1 = rep(NA, 20); t2 = rep(NA, 20)
t1[seq(1,20,by=2)] = 1:10
t2[seq(1,20,by=2)] = t1[seq(1,20,by=2)] + sample(1:5, 10, replace=T) + abs(rnorm(10))
t1[seq(2,20,by=2)] = t2[seq(1,20,by=2)]
t2[seq(2,20,by=2)] = t1[seq(2,20,by=2)] + sample(1:5) + abs(rnorm(10))

kitty = rep(1:10, e=2)
kittyblarg= t2 + rnorm(20, sd=5)
kittyhappy = rep(0:1, times=5, e=2)
die = 0:1
cens = c(0, 0)
kittydied = ifelse(runif(20)>=.5, die, cens)
d = data.frame(kitty, kittyblarg, kittyhappy, 
               t1, t2, kittydied)

# Inspect the Surv object if desired
# Surv(t1,t2, kittydied)

# Inspect the data
d


### Create a function to feed to optim ###
pl_tv = function(pars, preds, died, t1, t2, data) {
  # Same arguments as before though will take a data object
  # plus variable names via string input. Also requires beginning
  # and end time point (t1, t2)
  dat = data[,c(preds, died, t1, t2)]
  dat = dat[order(dat$t2),]
  b = pars
  X = as.matrix(dat[,preds])
  died2 = dat[,died]
  LP = X%*%b
  ll = numeric(nrow(X))
  rows = 1:nrow(dat)
  for (i in rows){
    st_i = dat$t2[i]
    riskset = ifelse(rows < i | dat$t1 > st_i, F, T)              # if they have already died/censored (row < i) or 
    ll[i] = died2[i]*(LP[i] - log(sum(exp(LP[riskset]))) )        # if the initial time is greater than current end time (t1 > st_i), 
  }                                                               # they are not in the risk set, else they are.
  -sum(ll)
}


### Run the model via optim and compare to survival ###
init.val = c(0, 0)
out = optim(par = init.val, fn = pl_tv, preds = c('kittyblarg','kittyhappy'),
            died = 'kittydied', data=d, t1='t1', t2='t2', 
            method="BFGS", hessian=T)
out

# create a summary table
B = out$par
se = sqrt(diag(solve(out$hessian)))
Z = B/se
data.frame(B, exp=exp(B), se, Z, 
           p=ifelse( Z>0, pnorm(Z, lower=F)*2, pnorm(Z, lower=T)*2) )

# Compare to survival package
coxmod2 = coxph(Surv(t1, t2, kittydied) ~ kittyblarg + kittyhappy, method='breslow', 
                control=coxph.control(iter.max=1000))
coxmod2; coxmod2$loglik[2]
# summary(coxmod2)



############################
### Stratified Cox Model ###
############################

### Get the data ###
attach(ovarian)


### Create a function to feed to optim ###
# requires pl function above though one could extend to pl_tv
pl_strat = function(pars, preds, died, t, strata) {
  strat = as.factor(strata)
  d = data.frame(preds, died, t, strat)
  dlist = split(d, strata)
  neglls = sapply(dlist, FUN = function(x) pl(pars=pars, preds=x[,colnames(preds)], 
                                              died=x$died, t=x$t))
  sum(neglls)
}


### Run the model via optim and compare to survival ###
init.val = c(0, 0)
out = optim(par=init.val, fn=pl_strat, preds=ovarian[,c('age','ecog.ps')],
            died=fustat, t=futime, strata=rx,
            method="BFGS", hessian=T)
out
B = out$par
se = sqrt(diag(solve(out$hessian)))
Z = B/se
data.frame(B, exp=exp(B), se, Z, 
           p=ifelse( Z>0, pnorm(Z, lower=F)*2, pnorm(Z, lower=T)*2) )

mod_compare = coxph(Surv(futime, fustat) ~ age + ecog.ps + strata(rx), data=ovarian)
mod_compare; mod_compare$loglik[2]
# summary(mod_compare)
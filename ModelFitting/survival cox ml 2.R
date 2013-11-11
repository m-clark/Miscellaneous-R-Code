#time var coef

### Create Data ###
set.seed(123)

t1 = rep(NA, 20); t2 = rep(NA, 20)
t1[seq(1,20,by=2)] = 1:10
t2[seq(1,20,by=2)] = t1[seq(1,20,by=2)]+sample(1:5,10, replace=T) + abs(rnorm(10))
t1[seq(2,20,by=2)] = t2[seq(1,20,by=2)]
t2[seq(2,20,by=2)] = t1[seq(2,20,by=2)]+sample(1:5) + abs(rnorm(10))

cbind(t1,t2)
Surv(t1,t2, kittydied)

kitty = rep(1:10, e=2)
kittyblarg= t2+rnorm(20, sd=5)
kittyhappy = rep(0:1,times=5, e=2)
die = 0:1
cens = c(0,0)
kittydied = ifelse(runif(20)>=.5, die, cens)
d = data.frame(kitty,kittyblarg,kittyhappy,t1,t2,kittydied)

d

#standard
init.val = c(0,0)
out = optim(par = init.val, fn = pl, preds = d[,c('kittyblarg','kittyhappy')],
            died = d[,'kittydied'], t=t2, #id = d[,'kitty'], 
             hessian=T) #method="Nelder",
out
B=out$par
se=diag(sqrt(solve(out$hessian)))
Z=B/se
data.frame(B, exp=exp(B),se, Z, p=ifelse(Z>=0,pnorm(Z, lower=F)*2,pnorm(Z)*2 ))

coxmod = coxph(Surv(t2, kittydied)~kittyblarg+kittyhappy, method='breslow', control=coxph.control(iter.max=1000))
coxmod; coxmod$loglik

#multientry
init.val = c(0,0)
out = optim(par = init.val, fn = pl2, preds = d[,c('kittyblarg','kittyhappy')],
            died = d[,'kittydied'], t1=t1,t2=t2, #id = d[,'kitty'], 
            hessian=T) #method="Nelder",
out
B=out$par
se=sqrt(diag(solve(out$hessian)))
Z=B/se
data.frame(B, exp=exp(B),se, Z, p=ifelse(Z>=0,pnorm(Z, lower=F)*2,pnorm(Z)*2 ))


coxmod2 = coxph(Surv(t1, t2, kittydied)~kittyblarg+kittyhappy, method='breslow', control=coxph.control(iter.max=1000))
coxmod2; coxmod2$loglik[2]
summary(coxmod2)


#multientry alternate
init.val = c(0,0)
out = optim(par = init.val, fn = pl2B, preds = c('kittyblarg','kittyhappy'),
            died = 'kittydied', data=d, t1='t1',t2='t2', #id = d[,'kitty'], 
            hessian=T) #method="Nelder",
out

library(survival)


attach(ovarian)

init.val = c(1,1)
out = optim(par = init.val, fn = pl_strat, preds = ovarian[,c('age','ecog.ps')],died = fustat, t=futime, strata=rx,
            method="L-BFGS-B", hessian=T)
out
B=out$par
se=sqrt(diag(solve(out$hessian)))
Z=B/se
data.frame(B, exp=exp(B),se, Z,p=ifelse(Z>=0,pnorm(Z, lower=F)*2,pnorm(Z)*2 ) )

mod_compare = coxph(Surv(futime, fustat) ~ age + ecog.ps + strata(rx), data=ovarian)
mod_compare; mod_compare$loglik[2]
summary(mod_compare)


# library(foreign)
# write.dta(d, "C:/Users/mclark19/Desktop/survtest.dta")


#standard
pl = function(pars, preds, died, t) {
  b = pars
  X = preds[order(t),]
  died2 = died[order(t)]
  LP = b%*%t(X)
  ll = vector()
  rows = 1:nrow(preds)
  for (i in rows){
    riskset = ifelse(rows < i ,F,T)
    ll[i] = died2[i]*(LP[i] - log(sum(exp(LP[riskset]))) )
  }
  -sum(ll)
}

#new func for interval notation
pl2 = function(pars, preds, died, t1, t2) {
  dat = cbind(preds,died,t1,t2)
  ord = order(t2)
  dat = dat[ord,]
  b = pars
  X = dat[,1:ncol(preds)]
  died2 = dat[,ncol(preds)+1]
  LP = b%*%t(X)
  ll = vector()
  rows = 1:nrow(dat)
  for (i in rows){
    st_i = dat$t2[i]
    riskset = ifelse(rows < i | dat$t1 > st_i,F,T)
    ll[i] = died2[i]*(LP[i] - log(sum(exp(LP[riskset]))) )
  }
  -sum(ll)
}


###  alternate 
pl2B = function(pars, preds, died, t1, t2, data) {
  dat = data[,c(preds,died,t1,t2)]
  dat = dat[order(dat$t2),]
  b = pars
  X = dat[,preds]
  died2 = dat[,died]
  LP = b%*%t(X)
  ll = matrix(NA, nrow=nrow(dat))
  rows = 1:nrow(dat)
  for (i in rows){
    st_i = dat$t2[i]
    riskset = ifelse(rows < i | dat$t1 > st_i,F,T)                #if they have already died/censored (row < i) or if the initial time is
    ll[i] = died2[i]*(LP[i] - log(sum(exp(LP[riskset]))) )        #greater than current end time (t1 > st_i), they are not in the risk set, else they are.
  }
  -sum(ll)
}






#stratified cox, requires pl1 but easily extended pl2
pl_strat = function(pars, preds, died, t, strata) {
  strat = as.factor(strata)
  dt = data.frame(preds, died,t,strat)
  dlist = split(dt,strata)
  lls = lapply(dlist,FUN= function(x) pl(pars=pars, preds=x[,colnames(preds)], died=x$died,t=x$t))
  sum(unlist(lls))
}


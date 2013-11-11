#survival likelihood

set.seed(12)

dur = 1:10
kittyblarg= rnorm(10)
kittyhappy = rep(0:1,times=5)
kittydied = sample(0:1,10,replace=T)
d = data.frame(kittyblarg,kittyhappy,dur,kittydied)

d

pl = function(pars, preds, died) {
  b = pars
  X = preds
  LP = b%*%t(X)
  ll = vector()
  rows = 1:nrow(preds)
  for (i in rows){
    riskset = ifelse(rows < i ,F,T)
    ll[i] = died[i]*(LP[i] - log(sum(exp(LP[riskset]))) )
  }
  -sum(ll)
}

init.val = c(1,1)
out = optim(par = init.val, fn = pl, preds = d[,c('kittyblarg','kittyhappy')],died = d[,'kittydied'], method="L-BFGS-B", hessian=T)
out
B = out$par
se = sqrt(diag(solve(out$hessian)))
Z = B/se
data.frame(B, exp=exp(B),se, Z, p=pnorm(Z, lower=F)*2 )
#exp(out$par)

library(survival)
coxmod = coxph(Surv(dur, kittydied)~kittyblarg+kittyhappy)
coxmod; coxmod$loglik[2]
summary(coxmod)


#time var coef

set.seed(12)

dur = rep(NA, 20)
dur[seq(1,20,by=2)] = 1:10
dur[seq(2,20,by=2)] = 2:11

kitty = rep(1:10, e=2)
kittyblarg = rnorm(20)
kittyhappy = rep(0:1,times=5, e=2)
die = 0:1
cens = c(0,0)
kittydied = ifelse(runif(20)>=.5, die, cens)
d = data.frame(kitty,kittyblarg,kittyhappy,dur,kittydied)

d

#new func
pl2 = pl

init.val = c(1,1)
out = optim(par = init.val, fn = pl2, preds = d[,c('kittyblarg','kittyhappy')],
            died = d[,'kittydied'], #id = d[,'kitty'], 
            method="L-BFGS-B", hessian=T)
out
B = out$par
se = sqrt(diag(solve(out$hessian)))
Z = B/se
data.frame(B, exp=exp(B),se, Z, p=pnorm(Z, lower=F)*2 )

library(survival)
coxmod = coxph(Surv(dur, kittydied)~kittyblarg+kittyhappy)
coxmod; coxmod$loglik[2]
summary(coxmod)




# library(foreign)
# write.dta(d, "C:/Users/mclark19/Desktop/survtest.dta")
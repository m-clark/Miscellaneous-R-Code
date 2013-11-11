library(mvtnorm)
library(psych)

#Create the data set
set.seed(123)

lambda = (matrix(c(1,.5,.3,.6,0,0,0,0,
                 0,0,0,0,1,.7,.4,.5), 
               nrow=2, byrow=T))

phi = matrix(c(1,.25,.25,1),nrow=2,byrow=T)  #correlation of factors

factors = rmvnorm(1000, mean=rep(0,2), sigma=phi, "chol")
e = rmvnorm(1000, sigma=diag(8))

y = 0 + factors%*%lambda + e

#dim(y)
describe(y)
round(cor(y),3)

fa(y,nfactors=2, rotate="oblimin") #example exploratory


#measurement model
cfa.cov <- function (parms, data) {
  S = cov(data)*((nrow(data)-1)/nrow(data))  #ML covariance div by N rather than N-1, the multiplier adjusts
  lambda = matrix(c(1,parms[1:3],rep(0,4),rep(0,4), 1,parms[4:6]), ncol=2)
  dist.init = parms[7:14]
  disturbs = diag(dist.init)
  phi.init = matrix(c(parms[15],parms[17],parms[17],parms[16]), ncol=2) 
  sigtheta = lambda%*%phi.init%*%t(lambda)+disturbs
  q = dim(data)[2]  #or in Bollen p + q (here just q) = tr(data)
  #out = -(log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - q )  #a reduced version; Bollen 1989 p.107
  out = ((-(nrow(data))*q/2)*log(2*pi)) - ((nrow(data))/2)* (log(det(sigtheta)) + tr(S%*%solve(sigtheta)))  #should be same as Mplus H0 loglike
  -out
}


cfa.cor = function (parms, data) {
  S = cor(data)
  lambda = matrix(c(parms[1:4],rep(0,4),rep(0,4), parms[5:8]), ncol=2)
  dist.init = parms[9:16]
  disturbs = diag(dist.init)
  phi.init = matrix(c(1,parms[17],parms[17],1), ncol=2)
  sigtheta = lambda%*%phi.init%*%t(lambda)+disturbs
  q = dim(data)[2]
  #out = (log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - q )
  out = ((-nrow(data)*q/2)*log(2*pi))-(nrow(data)/2)* (log(det(sigtheta)) + tr(S%*%solve(sigtheta)))
  -out
}


par.init.cor=c(rep(0,8), rep(.9,8), 0) #for cor
out.cor = optim(par = par.init.cor, fn=cfa.cor, data=y, method="L-BFGS-B", lower = 0, upper = 1, control = c(list(fnscale = 1))) #, control = c(list(fnscale = 1), parscale = rep(0.01, length(par.init)))
loads.cor= matrix(c(out.cor$par[1:4],rep(0,4), rep(0,4),out.cor$par[5:8]), ncol=2)
disturbs.cor = out.cor$par[9:16]
#round(data.frame(loads.cor,disturbs.cor),3); round(matrix(c(1,out.cor$par[c(17,17)],1), ncol=2),3)

par.init.cov=c(rep(0,6), rep(.9,8), rep(.5,3)) #cov
out.cov = optim(par = par.init.cov, fn=cfa.cov, data=y, method="L-BFGS-B", lower = 0, control = c(list(fnscale = 1))) #, control = c(list(fnscale = 1), parscale = rep(0.01, length(par.init)))
loads.cov= matrix(c(1,out.cov$par[1:3],rep(0,4), rep(0,4),1,out.cov$par[4:6]), ncol=2)
disturbs.cov = out.cov$par[7:14]
#round(data.frame(loads.cov,disturbs.cov),3); round(matrix(c(out.cov$par[c(15,17,17,16)]), ncol=2),3)

output=list(out.raw = round(data.frame(loads.cov,disturbs.cov),3), cov.fact = round(matrix(c(out.cov$par[c(15,17,17,16)]), ncol=2),3),
            out.std = round(data.frame(loads.cor,disturbs.cor,Rsq=(1-disturbs.cor)),3), cor.fact = round(matrix(c(1,out.cor$par[c(17,17)],1), ncol=2),3),
            ll = out.cov$value)
output

######################
### Confirm lavaan ###
######################

library(lavaan)
y = data.frame(y); colnames(y) <- paste('x',1:8, sep="")
model <- ' F1  =~ x1 + x2 + x3 + x4
           F2  =~ x5 + x6 + x7 + x8 '

fit <- cfa(model, data=y)  #std.lv=T for standardized
summary(fit, fit.measures=TRUE)



########################
### Confirm in Mplus ###
########################

library(MplusAutomation)
prepareMplusData(data.frame(y),"C:/Users/mclark19/Desktop/factsim.dat")

# MODEL:
#   F1 BY X1-X4;
#   F2 BY X5-X8;
# 
# OUTPUT:
#   STDYX;

# #measurement model test
# cfa = function (parms, data) {
#   S = cov(y)
#   lambda = matrix(c(par.init[1:4],rep(0,4),rep(0,4), par.init[5:8]), ncol=2)
#   #lambda = matrix(c(1,par.init[1:3],rep(0,4),rep(0,4), 1,par.init[4:6]), ncol=2)
#   dist.init = par.init[7:14]
#   disturbs = diag(dist.init)
#   phi.init = matrix(c(1,par.init[15],par.init[15],1), ncol=2)
#   #phi = lower.tri(phi.init, diag=T)
#   sigtheta = lambda%*%phi.init%*%t(lambda)+disturbs
#   q = dim(y)[2]
# 
#   out =  (log(det(sigtheta)) + tr(S%*%solve(sigtheta)) - log(det(S)) - q )
#   out
# }


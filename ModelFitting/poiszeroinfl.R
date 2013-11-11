#Zeroinf Pois
ZIP <- function(y, X, par) {
  logitpars = par[grep('logit', names(par))]
  poispars = par[grep('pois', names(par))]
  
  #Logit part; in the following Xlogit = Xpois but the code could is setup to have different predictor sets via an additional argument
  Xlogit <- X
  
  LPlogit <- Xlogit%*%logitpars
  logi0 <- plogis(LPlogit)  #alternative 1/(1+exp(-LPlogit))
    
  #Poisson part
  Xpois <- X
  ypois <- y
  mupois <- exp(Xpois%*%poispars)
  
  ##LLs; first based on pscl code
  logliklogit = log(logi0 + exp(log(1 - logi0) - mupois))
  loglikpois = log(1-logi0) + dpois(y, lambda = mupois, log = TRUE)
  
  #Hilbe formulation
#   logliklogit = log(logi0 + (1 - logi0)*exp(- mupois) )
#   loglikpois = log(1-logi0) -mupois + log(mupois)*y     #not necessary: - log(gamma(y+1))
    
  y0 = y==0  #0 values
  yc = y>0   #Count part

  loglik = -( sum(logliklogit[y0])+sum(loglikpois[yc]) )
  loglik
}


library(foreign); library(pscl)
fish <- read.dta("http://www.stata-press.com/data/r11/fish.dta")
init.mod <- glm(count~persons + livebait, data=fish, x=T,y=T, "poisson")

#for these functions, a named vector for the starting values; for zip: need logit, pois; for zinb: logit, negbin, theta;
starts = c(logit=coef(init.mod),pois=coef(init.mod))  

### Poisson zero-inflated ###
#note, use BFGS; have had notable issues with Nelder default
optPois1 <- optim(par=starts , fn=ZIP, X=init.mod$x, y=init.mod$y, method="BFGS",control=list(maxit=5000,reltol=1e-15), hessian=T)
optPois1

B = optPois1$par
se = sqrt(diag(solve((optPois1$hessian))))
Z = B/se
p= ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)

zipoismod <- zeroinfl(count~persons + livebait, data=fish, dist="poisson")  #control=list(maxit=5000),
summary(zipoismod)
data.frame(B,se, Z, p)



#Zeroinf NB
ZINBloglik <- function(y, X, par) {
  logitpars = par[grep('logit', names(par))]
  negbinpars = par[grep('negbin', names(par))]
  theta = exp(par[grep('theta', names(par))])
  
  #Logit part
  Xlogit <- X
  LPlogit <- Xlogit%*%logitpars
  logi0 <-  plogis(LPlogit) #alternative 1/(1+exp(-LPlogit))
  
  #Negbin part
  Xnegbin <- X
  munb <- exp(Xnegbin%*%negbinpars)
  
  #theta part for Hilbe formulation
  alpha = 1/theta  #this discerned via comparison of R/Stata output and boards
  m = 1/alpha
  p = 1/(1+ alpha*munb)
  
  ##LLs; first based on pscl code  
   logliklogit = log( logi0 + exp(log(1 - logi0) + suppressWarnings(dnbinom(0,size = theta, mu = munb, log = TRUE))) )
   logliknegbin = log(1 - logi0) + suppressWarnings(dnbinom(y, size = theta, mu = munb, log = TRUE))
  
  #Hilbe formulation, in progress; logi part ok
#   logliklogit = log( logi0 + (1 - logi0)*(p^m) )
#   logliknegbin = log(1-logi0) + log(gamma(m+y)) - log(gamma(m)) + m*log(p) + y*log(1-p)   # gamma(y+1) not needed
  
  y0 = y==0
  yc = y>0
  
  loglik = -( sum(logliklogit[y0])+sum(logliknegbin[yc]) )
  loglik
}



#for these functions, a named vector for the starting values; for zip: need logit, pois; for zinb: logit, negbin, theta;
init.mod <- glm(count~persons + livebait, data=fish, x=T,y=T, "poisson")
startlogi = glm((count==0)~persons + livebait, data=fish, family="binomial")
startcount = glm(count~persons + livebait, data=fish, family="quasipoisson")

starts = c(logit=coef(startlogi), negbin=coef(startcount), theta=1)  

### Negative Binomial zero-inflated ###
optNB1 <- optim(par=starts , fn=ZINBloglik, X=init.mod$x, y=init.mod$y, method="BFGS",control=list(maxit=5000,reltol=1e-15), hessian=T)
optNB1

B = optNB1$par
se=sqrt(diag(solve((optNB1$hessian))))
Z = B/se
p= ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)

zinbmod1 <- zeroinfl(count~persons + livebait, data=fish, dist="negbin")  #control=list(maxit=5000),
summary(zinbmod1)
data.frame(B,se, Z, p)  #note that theta here is actually log(theta)


data("bioChemists", package = "pscl")

#for these functions, a named vector for the starting values; for zip: need logit, pois; for zinb: logit, negbin, theta;
init.mod = glm(art~fem + mar + kid5 + phd + ment, data=bioChemists, family="quasipoisson", x=T, y=T)
startlogi = glm((art==0)~fem + mar + kid5 + phd + ment, data=bioChemists, family="binomial")
startcount = glm(art~fem + mar + kid5 + phd + ment, data=bioChemists, family="quasipoisson")

starts = c(logit=coef(startlogi), negbin=coef(startcount), theta=summary(startcount)$dispersion)  

### Negative Binomial zero-inflated ###
optNB2 <- optim(par=starts , fn=ZINBloglik, X=init.mod$x, y=init.mod$y, method="BFGS",control=list(maxit=5000,reltol=1e-15), hessian=T)
optNB2

B = optNB2$par
se=sqrt(diag(solve((optNB2$hessian))))
Z = B/se
p= ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)

zinbmod <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")  #control=list(maxit=5000),
summary(zinbmod)
data.frame(B,se, Z, p)



# library(foreign)
# write.dta(bioChemists, "C:/Users/mclark19/Desktop/miscdata/bioChemists.dta")
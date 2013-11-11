#run formula code at bottom first

### Negative Binomial zero-inflated ###
library(pscl)
data(bioChemists)
zinbmod <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin", x=T, y=T)  #control=list(maxit=5000),
summary(zinbmod)


#profile theta
starts = c(logit=zinbmod$coef$zero, negbin=zinbmod$coef$count)
logthetas = seq(.05,4,length=1000)
lls = ZINBloglik_profiletheta(par=starts ,X=zinbmod$x[[1]], y=zinbmod$y, thetalist=logthetas)
summary(lls)
library(ggplot2)
qplot(x=logthetas, y=lls, geom="line") +
  geom_abline(intercept = min(lls), slope = 0) +
  ylim(c(1540,1600)) + 
  ylab("Negative log likelihood") +
  annotate("text", x = .75, y = 1548, label=paste("Min = ",round(min(lls),2), sep="") )


#profile women kids

coefs = cbind(seq(-.5,-.00,length=250), seq(-.5,-.00,length=250))
starts = c(logit=zinbmod$coef$zero, negbin=zinbmod$coef$count, theta=0.9763565)
lls = ZINBloglik_profilekidwomen(par=starts ,X=zinbmod$x[[1]], y=zinbmod$y, coeflist=coefs)
lls = rbind(lls, c(-0.1955068313,-0.1517324581,1549.99))
colnames(lls) = c('woman','kids','negll')
summary(lls)
wlims = c(-0.1955037756+2*0.0755926 , -0.1955037756-2*0.0755926)  
klims = c(-0.1517324581 +2*0.0542061, -0.1517324581 -2*0.0542061)

square95 = c(rlc = c(wlims[1],klims[2]),
             ruc = c(wlims[1],klims[1]),
             llc = c(wlims[2],klims[2]),
             ruc = c(wlims[2],klims[1]))

square95 = data.frame(x=c(wlims[1],wlims[1],wlims[2],wlims[2]),
           y=c(klims[2],klims[1],klims[1],klims[2]))

library(ggplot2)
gdat = data.frame(lls)
ggplot(aes(x=woman, y=kids), data=lls) +
  geom_point(aes(color=negll, size=-negll) ) + #, position="jitter"
  scale_size_continuous(range=c(.01,2)) +
  scale_color_gradient(high='gray10',low="blue") +
  #geom_vline(xintercept=wlims) +
  #geom_hline(yintercept=klims) +
  geom_polygon(aes(x=x,y=y), fill="gray90", data=square95, alpha=.15) +
  geom_point(color="darkred", size=5, data=lls[nrow(lls),]) +
  theme_bw()


#### OPTIMATOR
library(optimx)
init.mod = glm(art~fem + mar + kid5 + phd + ment, data=bioChemists, family="quasipoisson", x=T, y=T)
startlogi = glm((art==0)~fem + mar + kid5 + phd + ment, data=bioChemists, family="binomial")
startcount = glm(art~fem + mar + kid5 + phd + ment, data=bioChemists, family="quasipoisson")

starts = c(logit=coef(startlogi), negbin=coef(startcount), theta=summary(startcount)$dispersion)
optxZINB <- optimx(par=starts , fn=ZINBloglik, X=init.mod$x, y=init.mod$y, method="BFGS", hessian=T,
                   control=list(maxit=5000,reltol=1e-15, save.failures=T))
str(optxZINB)
print(attr(optxZINB,"details"))

optxZINB$par

plot(optxZINB)

save.image("C:/Users/mclark19/Desktop/CSR/Clients/Me/Data Files/Working/zerofinflated_profilell.RData")

#Zeroinf NB
ZINBloglik_profiletheta <- function(y, X, par, thetalist) {
  loglik = vector()
  for (i in 1:length(thetalist)){
  theta = exp(thetalist[i])
  
  logitpars = par[grep('logit', names(par))]
  negbinpars = par[grep('negbin', names(par))]
  
  #Logit part
  Xlogit <- X
  LPlogit <- Xlogit%*%logitpars
  logi0 <-  plogis(LPlogit) #alternative 1/(1+exp(-LPlogit))
  
  #Negbin part
  Xnegbin <- X
  munb <- exp(Xnegbin%*%negbinpars)

  ##LLs; first based on pscl code  
  logliklogit = log( logi0 + exp(log(1 - logi0) + suppressWarnings(dnbinom(0,size = theta, mu = munb, log = TRUE))) )
  logliknegbin = log(1 - logi0) + suppressWarnings(dnbinom(y, size = theta, mu = munb, log = TRUE))
  
  y0 = y==0
  yc = y>0
  
  loglik[i] = -( sum(logliklogit[y0])+sum(logliknegbin[yc]) )
  }
  loglik
}




#Zeroinf NB
ZINBloglik_profilekidwomen <- function(y, X, par, coeflist) {
  loglik_mat = data.frame(expand.grid(coeflist[,1],coeflist[,2]),ll=NA)
  logitpars = par[grep('logit', names(par))]
  negbinpars.init = par[grep('negbin', names(par))]
  theta = exp(par[grep('theta', names(par))])
  
  llseq = 1:nrow(loglik_mat)
  
    for(k in llseq){

        woman = loglik_mat[k,1]
        kids = loglik_mat[k,2]
        negbinpars = c( negbinpars.init[1], woman, negbinpars.init[3], kids, negbinpars.init[5:6])
        
        #Logit part
        Xlogit <- X
        LPlogit <- Xlogit%*%logitpars
        logi0 <-  plogis(LPlogit) #alternative 1/(1+exp(-LPlogit))
        
        #Negbin part
        Xnegbin <- X
        munb <- exp(Xnegbin%*%negbinpars)
        
        ##LLs; first based on pscl code  
        logliklogit = log( logi0 + exp(log(1 - logi0) + suppressWarnings(dnbinom(0,size = theta, mu = munb, log = TRUE))) )
        logliknegbin = log(1 - logi0) + suppressWarnings(dnbinom(y, size = theta, mu = munb, log = TRUE))
        
        y0 = y==0
        yc = y>0
      
        loglik_mat[k,"ll"] = -( sum(logliklogit[y0])+sum(logliknegbin[yc]) )
    }
  loglik_mat
}
       
       
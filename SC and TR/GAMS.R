#-------------------------------------------------------------------------------#
# See http://www.nd.edu/~mclark19/learn/GAMS.pdf for the handout on generalized #
# additive models this code regards.                                            #
#-------------------------------------------------------------------------------#

## mymod = lm(y ~ x1 + x2, data=mydata)



d = read.csv('http://www.nd.edu/~mclark19/learn/data/pisasci2006.csv')
library(psych)
describe(d)[-1,1:9]  #univariate



## library(car)
## scatterplotMatrix(d[,-c(1,3:5)],pch=19,cex=.5,reg.line=F, lwd.smooth=1.25,
##                   spread=F,ellipse=T, col=c('gray60','#2957FF','#FF8000'),
##                   col.axis='gray50')



## library(ggplot2); library(reshape2)
## #get data into a form to take advantage of ggplot
## dmelt = melt(d, id=c('Country','Overall'),
##              measure=c('Interest','Support','Income','Health','Edu','HDI'))
## 
## #leave the smooth off for now
## ggplot(aes(x=value,y=Overall), data=dmelt) +
##   geom_point(color='#FF8000',alpha=.75) +
##   #geom_smooth(se=F) +
##   geom_text(aes(label=Country), alpha=.25, size=1,angle=30, hjust=-.2,
##             vjust=-.2) +
##   facet_wrap(~variable, scales='free_x') +
##   ggtheme



## library(mgcv)
## ggplot(aes(x=value,y=Overall), data=dmelt) +
##   geom_point(color='#FF8000',alpha=.75) +
##   geom_smooth(se=F, method='gam', formula=y~s(x), color='#2957FF') +
##   facet_wrap(~variable, scales='free_x') +
##   ggtheme



library(mgcv)
mod_lm <- gam(Overall ~ Income, data=d)
summary(mod_lm)



mod_gam1 <- gam(Overall ~ s(Income, bs="cr"), data=d)
summary(mod_gam1)



## plot(mod_gam1)



AIC(mod_lm)
summary(mod_lm)$sp.criterion
summary(mod_lm)$r.sq  #adjusted R squared



anova(mod_lm, mod_gam1, test="Chisq")



mod_lm2 <- gam(Overall ~ Income + Edu + Health, data=d)
summary(mod_lm2)



mod_gam2 <- gam(Overall ~ s(Income) + s(Edu) + s(Health), data=d)
summary(mod_gam2)



mod_gam2B = update(mod_gam2, .~.-s(Health) + Health)
summary(mod_gam2B)



## plot(mod_gam2, pages=1, residuals=T, pch=19, cex=0.25,
##      scheme=1, col='#FF8000', shade=T,shade.col='gray90')



## # Note that mod_gam2$model is the data that was used in the modeling process,
## # so it will have NAs removed.
## testdata = data.frame(Income=seq(.4,1, length=100),
##                       Edu=mean(mod_gam2$model$Edu),
##                       Health=mean(mod_gam2$model$Health))
## fits = predict(mod_gam2, newdata=testdata, type='response', se=T)
## predicts = data.frame(testdata, fits)
## 
## ggplot(aes(x=Income,y=fit), data=predicts) +
##   geom_smooth(aes(ymin = fit - 1.96*se.fit, ymax=fit + 1.96*se.fit),
##               fill='gray80', size=1,stat='identity') +
##   ggtheme



## vis.gam(mod_gam2, type='response', plot.type='contour')



mod_gam3 <- gam(Overall ~ te(Income,Edu), data=d)
summary(mod_gam3)



## vis.gam(mod_gam3, type='response', plot.type='persp',
##         phi=30, theta=30,n.grid=500, border=NA)



anova(mod_lm2, mod_gam2, test="Chisq")



gam.check(mod_gam2, k.rep=1000)



mod_1d = gam(Overall ~ s(Income) + s(Edu), data=d)
mod_2d = gam(Overall ~ te(Income,Edu, bs="tp"), data=d)
AIC(mod_1d, mod_2d)



mod_A = gam(Overall ~ s(Income, bs="cr", k=5) + s(Edu, bs="cr", k=5), data=d)
mod_B = gam(Overall ~ s(Income, bs="cr", k=5) + s(Edu, bs="cr", k=5) + te(Income, Edu), data=d)

anova(mod_A,mod_B, test="Chi")



## #set default options in case we use ggplot later
## ggtheme =
##   theme(
##     axis.text.x = element_text(colour='gray50'),
##     axis.text.y = element_text(colour='gray50'),
##     panel.background = element_blank(),
##     panel.grid.minor = element_blank(),
##     panel.grid.major = element_blank(),
##     panel.border = element_rect(colour='gray50'),
##     strip.background = element_blank()
##   )



## ggplot(aes(x=Income, y=Overall), data=d) +
##   geom_point(color="#FF8000") +
##   geom_smooth(se=F, method='gam', formula=y~s(x, bs="cr")) +
##   xlim(.4,1) +
##   ggtheme



## ############################
## ### Wood by-hand example ###
## ############################
## 
## 
## size = c(1.42,1.58,1.78,1.99,1.99,1.99,2.13,2.13,2.13,
##          2.32,2.32,2.32,2.32,2.32,2.43,2.43,2.78,2.98,2.98)
## wear = c(4.0,4.2,2.5,2.6,2.8,2.4,3.2,2.4,2.6,4.8,2.9,
##          3.8,3.0,2.7,3.1,3.3,3.0,2.8,1.7)
## x= size-min(size); x = x/max(x)
## d = data.frame(wear, x)
## 
## #cubic spline function
## rk <- function(x,z) {
##   ((z-0.5)^2 - 1/12)*((x-0.5)^2 - 1/12)/4-
##     ((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2 + 7/240) / 24
## }
## 
## spl.X <- function(x,knots){
##   q <- length(knots) + 2  # number of parameters
##   n <- length(x)          # number of observations
##   X <- matrix(1,n,q)      # initialized model matrix
##   X[,2] <- x              # set second column to x
##   X[,3:q] <- outer(x,knots,FUN=rk) # remaining to cubic spline
##   X
## }
## 
## spl.S <- function(knots) {
##   q = length(knots) + 2
##   S = matrix(0,q,q)     # initialize matrix
##   S[3:q,3:q] = outer(knots,knots,FUN=rk) # fill in non-zero part
##   S
## }
## 
## #matrix square root function
## mat.sqrt <- function(S){
##   d = eigen(S, symmetric=T)
##   rS = d$vectors%*%diag(d$values^.5)%*%t(d$vectors)
##   rS
## }
## 
## #the fitting function
## prs.fit <- function(y,x,knots,lambda){
##   q = length(knots) + 2   # dimension of basis
##   n = length(x)           # number of observations
##   Xa = rbind(spl.X(x,knots), mat.sqrt(spl.S(knots))*sqrt(lambda))  # augmented model matrix
##   y[(n+1):(n+q)] = 0   #augment the data vector
##   lm(y ~ Xa-1)   # fit and return penalized regression spline
## }



## knots = 1:4/5
## X = spl.X(x,knots)     # generate model matrix
## mod.1 = lm(wear~X-1)   # fit model
## xp <- 0:100/100        # x values for prediction
## Xp <- spl.X(xp, knots) # prediction matrix
## 
## # Base R plot
## plot(x, wear, xlab='Scaled Engine size', ylab='Wear Index', pch=19,
##      col="#FF8000", cex=.75, col.axis='gray50')
## lines(xp,Xp%*%coef(mod.1), col='#2957FF')   #plot
## 
## # ggplot
## library(ggplot2)
## ggplot(aes(x=x, y=wear), data=data.frame(x,wear))+
##   geom_point(color="#FF8000") +
##   geom_line(aes(x=xp, y=Xp%*%coef(mod.1)), data=data.frame(xp,Xp), color="#2957FF") +
##   ggtheme
## 



## knots = 1:7/8
## 
## d2 = data.frame(x=xp)
## 
## for (i in c(.1,.01,.001,.0001,.00001,.000001)){
##   mod.2 = prs.fit(wear, x, knots,i)   #fit penalized regression
##                                       #spline choosing lambda
##   Xp = spl.X(xp,knots) #matrix to map parameters to fitted values at xp
##   d2[,paste('lambda = ',i, sep="")] = Xp%*%coef(mod.2)
## }
## 
## ### ggplot
## library(ggplot2); library(reshape)
## d3 = melt(d2, id='x')
## 
## ggplot(aes(x=x, y=wear), data=d) +
##   geom_point(col='#FF8000') +
##   geom_line(aes(x=x,y=value), col="#2957FF", data=d3) +
##   facet_wrap(~variable) +
##   ggtheme
## 
## ### Base R approach
## par(mfrow=c(2,3))
## for (i in c(.1,.01,.001,.0001,.00001,.000001)){
##   mod.2 = prs.fit(wear, x, knots,i)
##   Xp = spl.X(xp,knots)
##   plot(x,wear, main=paste('lambda = ',i), pch=19,
##        col="#FF8000", cex=.75, col.axis='gray50')
##   lines(xp,Xp%*%coef(mod.2), col='#2957FF')
## }



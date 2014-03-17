#----------------------------------------------------------------------------------#
# The following provides a simple working example of a standard regression model   #
# using BUGS/R2OpenBugs to comppare with Stan/rStan. It is just for demonstration, # 
# and hopefully to allow some to #more easily jump right in to Bayesian methods    #
# if they are comfortable with R.  See rstan_linreg for comparison to Stan.        #
#----------------------------------------------------------------------------------#


#######################
### Create the Data ###
#######################

### create a correlation matrix of one's choosing assuming response as last column/row ###

cormat = matrix(c(1, .2, -.1, .3,
                  .2, 1, .1, .2,
                  -.1, .1, 1, .1,
                  .3, .2, .1, 1),
                ncol=4, byrow=T)

cormat

# ensure pos def
library(Matrix)
cormat = nearPD(cormat, corr=T)$mat

### generate data ###

library(MASS)
means = rep(0, ncol(cormat))
N = 1000
d = mvrnorm(N, means, cormat, empirical=T)
colnames(d) = c('X1', 'X2', 'X3', 'y')
d[,'y'] = d[,'y'] -.1 # unnecessary, just to model a non-zero intercept
str(d)
cor(d)


### prepare for later processing ###

# strip X (add intercept column) and y
X = cbind(1, d[,1:3]); colnames(X) = c('Intercept', 'X1', 'X2', 'X3')
y = d[,4]
K = ncol(X)

# for comparison
modlm = lm(y~., data.frame(d))



##################
### BUGS setup ###
##################

bugsdat = list('y', 'X', 'N', 'K')

sink('lmbugs.txt')
cat(
'model {
  for (n in 1:N){
    mu[n] <- beta[1]*X[n,1] + beta[2]*X[n,2] + beta[3]*X[n,3] + beta[4]*X[n,4]
    y[n] ~ dnorm(mu[n], tau.y)
  }
  
  for (k in 1:K){
    beta[k] ~ dnorm(0, .001)                                                    # prior for reg coefs
  }
  
  # Attempt at half-cauchy
  # Scale parameter is 5, so precision of xiNs = 1/5^2 = 0.04
  xi ~ dnorm(0, .04)I(0.001,)
  chSq ~ dgamma(0.5, 0.5)                                                       # chi^2 with 1 d.f.
  sigma.y <- xi/sqrt(chSq)                                                      # prior for sigma; cauchy = normal/sqrt(chi^2)
  tau.y <- pow(sigma.y, -2)                                                     # precision
}'
)
sink()


# initial values not necessary, but can specify as follows
# inits <- function(){
#   list(beta=rep(0,4), sigma.y=runif(1,0,10), xi=rnorm(1), tau.eta=runif(1) )
# }
parameters <- c('beta', 'sigma.y')



#####################
### Run the model ###
#####################
library(R2OpenBUGS)

# note that n.thin argument is used differently than other packages.  The 
# gist is, specify the n posterior draws (per chain) you to keep want as n.iter-n.burnin.  
# The thinned samples aren't stored. In other packages: n.iter is total before 
# thinning and including burnin, and n.keep is (n.iter-n.burnin) / n.thin.  With bugs, 
# n.keep is the same, but as far as arguments your you'll want to think of n.iter 
# your n posterior draws after thinning. So the following all produce 1000 posterior  
# draws in R2OpenBUGS:
# n.iter=3000, n.thin=1, n.burnin=2000
# n.iter=3000, n.thin=10, n.burnin=2000
# n.iter=3000, n.thin=100, n.burnin=2000
# 
# In other packages, with those arguments you'd end up with 1000, 300, 30 n posterior
# draws.

lmbugs <- bugs(bugsdat, inits=NULL, parameters, model.file='lmbugs.txt', n.chains=3, 
               n.iter=3000, n.thin=10, n.burnin=2000, codaPkg=F, debug=F)
print(lmbugs, digits=3)
plot(lmbugs)

library(coda); library(scales); library(ggthemes)
lmbugscoda = as.mcmc.list(lmbugs)
traceplot(lmbugscoda, col=alpha(gg_color_hue(3), .5))
densityplot(lmbugscoda, col=alpha(gg_color_hue(3), .5))
plot(lmbugscoda, col=alpha(gg_color_hue(3), .25))
corrplot:::corrplot(cor(lmbugscoda[[2]]))  # noticeably better than levelplot



############################
### Other visualizations ###
############################

### Playing with denstrip for other visuals
library(denstrip); library(scales)

betas = lmbugs$sims.list$beta
betameans = lmbugs$mean$beta
betasds = lmbugs$sd$beta
dens = sapply(1:4, function(i) dnorm(betas[,i], betameans[i], betasds[i]))

# initialize plot
plot(betas, xlim=c(-.3, .4), ylim=c(0, 5), xlab='Coefficients',  type='n', bty='n', ylab='', yaxt='n',
     col.lab='gray25')
axis(side=1, col='gray50', col.ticks='gray75', col.axis='gray25')
# using density from above
sapply(1:4, function(i) denstrip(betas[,i], at=i, width=.5, colmax=gg_color_hue(ncol(betas))[i]))

# kernel density
sapply(1:4, function(i) denstrip(betas[,i], dens[,i], at=i, width=.5))

# varying width
plot(betas, xlim=c(-.3, .4), ylim=c(0, 5), xlab='Coefficients',  type='n', bty='n', ylab='', yaxt='n',
     col.lab='gray25')
axis(side=1, col='gray50', col.ticks='gray75', col.axis='gray25')
sapply(1:4, function(i) vwstrip(betas[,i], dens[,i], at=i, width=.5, 
                                col=alpha(gg_color_hue(ncol(betas))[i], .5),
                                border=NA))


#-----------------------------------------------------------------------------------#
# The following provides a simple working example of a standard regression model    #
# using JAGS/rjags to comppare with Stan/rStan. It is just for demonstration,       # 
# and hopefully to allow some to more easily jump right in to Bayesian methods      #
# if they are comfortable with R.  See rstan_linreg and bugs_linreg for comparison. #
#-----------------------------------------------------------------------------------#


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

jagsdat = list('y'=y, 'X'=X, 'N'=N, 'K'=K)

# Note that for this sink part and the same section on the corresponding bugs code,
# you must run it intently. In other words, you cannot just ctrl+shft+enter the 
# whole script. See also 'ggsave'.
sink('lmjags.txt')
cat(
'model {
  for (n in 1:N){
    mu[n] <- beta[1]*X[n,1] + beta[2]*X[n,2] + beta[3]*X[n,3] + beta[4]*X[n,4]
    y[n] ~ dnorm(mu[n], inv.sigma.sq)
  }
  
  for (k in 1:K){
    beta[k] ~ dnorm(0, .001)                                                    # prior for reg coefs
  }
  
  # Half-cauchy as in Gelman 2006
  # Scale parameter is 5, so precision of z = 1/5^2 = 0.04
  z ~ dnorm(0, .04)I(0,)
  chSq ~ dgamma(0.5, 0.5)                                                       # chi^2 with 1 d.f.
  sigma.y <- z/sqrt(chSq)                                                       # prior for sigma; cauchy = normal/sqrt(chi^2)
  inv.sigma.sq <- pow(sigma.y, -2)                                              # precision
}'
)
sink()


# explicitly provided initial values not necessary, but can specify as follows
# inits <- function(){
#   list(beta=rep(0,4), sigma.y=runif(1,0,10) )
# }
parameters <- c('beta', 'sigma.y')



#####################
### Run the model ###
#####################
library(rjags)
lmjagsmod <- jags.model(file='lmjags.txt', data=jagsdat, # inits=inits
                        n.chains=3, n.adapt=2000)

# update(lmjagsmod, 10000)

lmjags = coda.samples(lmjagsmod, c('beta', 'sigma.y'), n.iter=10000, thin=10, n.chains=3)
summary(lmjags)
effectiveSize(lmjags)

# visulize
library(coda); library(scales); library(ggthemes)
traceplot(lmjags, col=alpha(gg_color_hue(3), .5))
densityplot(lmjags, col=alpha(gg_color_hue(3), .5))
plot(lmjags, col=alpha(gg_color_hue(3), .25))
corrplot:::corrplot(cor(lmjags[[2]]))  # noticeably better than levelplot
par(mar=c(5, 4, 4, 2) + 0.1) # reset margins


############################
### Other visualizations ###
############################

### Playing with denstrip for other visuals
library(denstrip); library(scales)

betas = do.call(rbind, lmjags)[,1:4]
betameans = summary(lmjags)$statistics[1:4,'Mean']
betasds = summary(lmjags)$statistics[1:4, 'SD']
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


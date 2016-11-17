#------------------------------------------------------------------------------------------#
# The following provides an example of beta regression using Stan/rstan, with comparison   #
# to results with R's betareg package.  Several data sets from that package are available, #
# to play with, but as they are a bit problematic in one way or another I also provide     #
# a simple simulated data set.  One can find another example or two of this type of model  #
# on the stan user boards.                                                                 #
#------------------------------------------------------------------------------------------#


#######################
### Create the Data ###
#######################

library(betareg)
# Data for assessing the contribution of non-verbal IQ to children's reading skills in dyslexic and non-dyslexic children.
# issue: 30% of data has a value of .99
# data("ReadingSkills")
# ?ReadingSkills
# y = ReadingSkills$accuracy 
# 
# brmod = betareg(accuracy ~ dyslexia + iq, data = ReadingSkills)
# X = cbind(1, scale(model.matrix(brmod)[,c('dyslexia','iq')], scale=F))


# or this, issue: ignores batch effects
# data("GasolineYield")
# ?GasolineYield
# y = GasolineYield$yield
# X = cbind(1, scale(GasolineYield[,c('gravity','pressure','temp')]))

# yet another data option, issue: only two binary predictors
# data(WeatherTask)
# ?WeatherTask
# y = WeatherTask$agreement
# brmod = betareg(agreement ~ priming + eliciting, data = WeatherTask)
# X = model.matrix(brmod)

# simulated data; probably a better illustration, or at least better behaved one.
set.seed(1234)

N = 500
x1 = rnorm(N)
x2 = rnorm(N)
X = cbind(1, x1, x2)
beta = c(-1,.2,-.3)
mu = plogis(X%*%beta)  # add noise if desired + rnorm(N, sd=.01)
phi = 10
A = mu*phi
B = (1-mu)*phi
y = rbeta(N, A, B)
hist(y, 'FD')

# model for later comparison
brmod = betareg(y ~ ., data = data.frame(y, X[,-1]))
summary(brmod)


# Stan data list
dat = list(N = length(y), K=ncol(X), y=y, X=X)



##################
### Stan Model ###
##################
stanmodelcode = '
data {
  int<lower=1> N;                   // sample size
  int<lower=1> K;                   // K predictors
  vector<lower=0,upper=1>[N] y;     // response 
  matrix[N,K] X;                    // predictor matrix
}

parameters {
  vector[K] theta;                  // reg coefficients
  real<lower=0> phi;                // dispersion parameter
}

transformed parameters{
  vector[K] beta;

  beta = theta * 10;               // "matt trick" if desired
}

model {
  // model calculations
  vector[N] Xbeta;                  // linear predictor
  vector[N] mu;                     // transformed linear predictor
  vector[N] A;             // parameter for beta distn
  vector[N] B;             // parameter for beta distn

  Xbeta = X * beta;
  for (i in 1:N) { 
    mu[i] = inv_logit(Xbeta[i]);   
  }

  A = mu * phi;
  B = (1.0 - mu) * phi;

  // priors
  theta ~ normal(0, 1);   
  phi ~ cauchy(0, 5);               // different options for phi  
  //phi ~ inv_gamma(.001, .001);
  //phi ~ uniform(0, 500);          // put upper on phi if using this

  // likelihood
  y ~ beta(A, B);
}
'

################
### Test Run ###
################
# Note you will get some informational messages but they appear to be inconsequential;

library(rstan)
test = stan(model_code = stanmodelcode, model_name = "betareg", #init=0, init_r=.01,
             pars = c('beta','phi'), data = dat, iter = 2000, 
             warmup=200, thin=1, chains = 3, verbose = F) 
print(test, digits=3)
# traceplot(test, inc_warmup=T)
# traceplot(test, inc_warmup=F)

summary(brmod)


#################
### Final Run ###
#################
fit = stan(model_code = stanmodelcode, model_name = "betareg", fit=test,
            data = dat, iter = 22000, warmup=2000, thin=20, 
            cores = 4, verbose = F) 
print(fit, pars = c('beta','phi'),digits_summary=3, probs = c(0, .025, .5, .975, 1))

summary(brmod)


# diagnostics
shinystan::launch_shinystan(fit)


### posterior predictive check ###
# extract posterior draws
betareg_sim = extract(fit, permuted=T)
nSim = length(betareg_sim$lp__)
yRep = array(NA, c(nSim, dat$N))

# simulate yRep as in BDA appendix C; could also do in generated quantities block of stan code
for (s in 1:nSim){
  betaS = betareg_sim$beta[s,]
  muS = plogis(X%*%betaS)
  A = muS*betareg_sim$phi[s]
  B = (1-muS)*betareg_sim$phi[s]
  yRep[s,] = rbeta(dat$N, A, B)
}
str(yRep)


# example, as in BDA p.144, code from appendix C; probably not very useful
par(mfrow=c(5,4), mar=c(4,4,2,2))
hist(y, xlab='', main='y', 'FD')
for (s in 1:19) hist(yRep[s,], xlab='', main=paste('yRep',s), 'FD')
layout(1)

# the following plots get at the same thing
library(ggplot2); library(reshape2)
gdat = data.frame(X,y)
gdatyRep = melt(yRep)

# individual observation distibutions; note the theme is my own (available on github), but otherwise you'll have to define yours.
ggplot() + 
  geom_line(aes(x=value, group=factor(Var2)), color='darkred', stat='density', lwd=.5, alpha=.05, data=gdatyRep) +
  stat_density(aes(x=y), data=gdat, alpha=.15, geom='density', color=NA, fill='black') +
  lazerhawk::theme_trueMinimal()

# whole response distributions, comparison with betareg fits (takes awhile)
ggplot(aes(x=y), xlim=c(.2,1), data=gdat) + 
  geom_line(aes(x=value,  group=as.factor(Var1), col=as.factor(Var1)), stat='density', data=gdatyRep, show_guide=F, alpha=.05) +
  stat_density(geom='line', lwd=2, alpha=.2) + 
  stat_density(aes(x=brmod$fitted), geom='line', col='darkred', lwd=2, alpha=.2) +
  stat_density(aes(x=value), geom='line', col='#FF5500', lwd=2, alpha=.2, data=gdatyRep) +
  lazerhawk::theme_trueMinimal()

# Examine quantiles
sapply(list(yRep=gdatyRep$value, brFitted=brmod$fitted.values, y=y), quantile, p=c(0,.1,.25,.5,.625,.75,.8,.9,.95,1))

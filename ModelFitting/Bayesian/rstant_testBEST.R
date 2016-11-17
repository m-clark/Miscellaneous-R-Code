#-------------------------------------------------------------------------------#
# The following is based on Kruschke's 2012 JEP article 'Bayesian estimation    #
# supercedes the t-test (BEST)' with only minor changes to stan model. It uses  #
# the JAGS/BUGS code in the paper's Appendix B as the reference.                #
#-------------------------------------------------------------------------------#

#######################
### Create the Data ###
#######################

### play around with the specs if you like
set.seed(1234)
N1 = 50
N2 = 50
mu1 = 1
mu2 = -.5
sig1 = 1
sig2 = 1
Ng = 2

y1 = rnorm(N1, mu1, sig1)
y2 = rnorm(N2, mu2, sig2)
y = c(y1, y2)

groupID = as.numeric(gl(2, N1))

## if unbalanced
# group = 1:2
# groupID = rep(group, c(N1,N2))

tapply(y, groupID, psych:::describe)

##################
### Stan Setup ###
##################
standat= list(N=length(y), Ng=Ng, groupID=groupID, y=y)

stanmodelcode = '
data {
  int<lower=1> N;                               // sample size (note:putting bounds provides simple data check)
  int<lower=2> Ng;                              // number of groups
  vector[N] y;                                  // response
  int<lower=1, upper=Ng> groupID[N];            // group ID
}

transformed data{
  real meany;                                   // mean of y; see mu prior
  
  meany = mean(y); 
}

parameters {
  vector[2] mu;                                 // estimated group means and sd
  vector<lower=0>[2] sigma;                     // Kruschke puts upper bound as well; ignored here
  real<lower=0, upper=100> nu;                  // df for t distribution
}

transformed parameters {                        // none needed

}

model {
  // priors
  mu ~ normal(meany, 10);                       // note that there is a faster implementation of this for stan; sd here is more informative than in Kruschke paper
  sigma ~ cauchy(0, 5);
  nu ~ exponential(1.0/29);                     // Based on Kruschke; makes mean nu 29 (might consider upper bound, too large and might as well switch to normal)
  
  
  // likelihood
  for (n in 1:N){
    y[n] ~ student_t(nu, mu[groupID[n]], sigma[groupID[n]]);
    //y[n] ~ normal(mu[groupID[n]], sigma[groupID[n]]);           // for comparison, remove all nu specifications if you do this
  }
}

generated quantities {
  vector[N] yRep;                               // posterior predictive distribution
  real muDiff;                                  // mean difference
  real CohensD;                                 // effect size; see footnote 1 in Kruschke paper
  real CLES;                                    // common language effect size
  real CLES2;                                   // a more explicit approach; the mean should roughly equal CLES

  for (n in 1:N){
    yRep[n] = student_t_rng(nu, mu[groupID[n]], sigma[groupID[n]]);
  }

  muDiff = mu[1] - mu[2];
  CohensD = muDiff / sqrt(sum(sigma)/2);
  CLES = normal_cdf(muDiff / sqrt(sum(sigma)), 0, 1);
  CLES2 = student_t_rng(nu, mu[1], sigma[1]) - student_t_rng(nu, mu[2], sigma[2]) > 0;
}

'



#############################
### Run and inspect model ###
#############################

### Run model/examine basic diagnostic plots
library(rstan) 
# you can ignore the informational message
fit = stan(model_code=stanmodelcode, data=standat, iter=12000, warmup=2000, cores=4, thin=10)

shinystan::launch_shinystan(fit)

### Print summary of model
print(fit, digits=3, pars=c('mu', 'sigma', 'muDiff', 'CohensD', 'CLES', 'CLES2','nu','lp__'))


### Extract quantities of interest for more processing/visualization.
yRep = extract(fit, par='yRep')$yRep

# compare population and observed data values to estimates in summary print
# mean difference
muDiff = extract(fit, par='muDiff')$muDiff
means = tapply(y, groupID, mean)
sds = tapply(y, groupID, sd)
mu1-mu2           # based on population values
abs(diff(means))  # observed in data

# Cohen's d
CohensD = extract(fit, par='CohensD')$CohensD
(mu1-mu2) / sqrt((sig1^2+sig2^2)/2)      # population
(means[1]-means[2]) / sqrt(sum(sds^2)/2)   # observed

# common language effect size is the probability that a randomly selected score from one 
# population will be greater than a randomly sampled score from the other
CLES = extract(fit, par='CLES')$CLES
pnorm((mu1-mu2) / sqrt(sig1^2+sig2^2))        # population
pnorm((means[1]-means[2]) / sqrt(sum(sds^2))) # observed



########################
### Model Comparison ###
########################

### Compare to Welch's t-test
t.test(y1,y2)


### Compare to BEST; note that it requires coda, whose traceplot function will overwrite rstan's
library(BEST)
BESTout = BESTmcmc(y1, y2, numSavedSteps=12000, thinSteps=10, burnInSteps=2000) 
summary(BESTout)



#####################
### Visualization ###
#####################

library(ggplot2); library(reshape2);

### plot posterior predictive distribution vs. observed data density
gdat = melt(yRep)
str(gdat)
colnames(gdat) = c('iteration', 'observation', 'value' )
gdat$groupID = factor(rep(groupID, e=2000))  # change this to match your sample size/chain length
gdat$observation = factor(gdat$observation)

ggplot(aes(x=value), data=gdat) + 
  geom_density(aes(group=groupID, fill=groupID), color=NA, alpha=.25) +
  geom_line(aes(group=observation, color=groupID), stat='density', alpha=.05) +
  geom_point(aes(x=y, y=0, color=factor(groupID)), alpha=.15, size=5, data=data.frame(y, groupID)) +
  xlim(c(-8,8)) +   # might get a warning if extreme values are cut out
  geom_density(aes(group=groupID, color=groupID, x=y), alpha=.05, data.frame(groupID=factor(groupID),y)) 

### plot mean difference or other values of interest
ggplot(aes(x=muDiff), data=data.frame(muDiff=muDiff)) + 
  geom_density(alpha=.25) +
  xlim(c(0,3.5)) +
  geom_point(x=muDiff, y=0, alpha=.01, size=3) +
  geom_path(aes(x=quantile(muDiff, c(.025, .975)), y=c(.2,.2)), size=2, alpha=.5, color='darkred', data=data.frame())


### BEST plots
par(mfrow=c(2,2))
sapply(c("mean", "sd", "effect", "nu"), function(p) plot(BESTout, which=p))
layout(1)

### Mixture model (draft); The Stan manual notes issues with mixture modeling, esp. label switching 
### (as does the bugs book),  but the code in the manual chapter on mixture 
### modeling does not appear able to reproduce the old faithful mixture or even 
### a simulated mixture of normals, even with some heavy hand holding on 
# starting points, limits etc.  On occasion the means can be recovered, but the 
# variances are all over the place and notably influenced by the prior.

data(faithful)
head(faithful)
par(mfrow=c(1,2))
apply(faithful, 2, hist, 'FD', freq=F)
layout(1)

y1 = rnorm(500, 50, 5)
y2 = rnorm(500, 80, 10)
y = c(y1, y2)
psych::describe(data.frame(y1, y2))

# take your pick
standat = list(N=nrow(faithful), K=2, y=faithful$waiting,
               mustart=quantile(faithful$waiting, prob=c(.25,.75)),
               sigmaupper=10, sigmastart=1)

standat = list(N=nrow(faithful), K=2, y=faithful$eruptions, 
               mustart=quantile(faithful$eruptions, prob=c(.25,.75)),
               sigmaupper=1, sigmastart=1)

standat = list(N=length(y), K=2, y=y, 
               mustart=quantile(y, prob=c(.1,.9)),       
               sigmaupper=sd(y), sigmastart=5)


stanmodelcode = '
data {
  int<lower=1> N;                             // Sample size
  int<lower=2> K;                             // K mixture components
  vector[N] y;                                // response
  vector[K] mustart;                          // starting values
  real sigmastart;                            // sigma to go with starting values
  real sigmaupper;                            // upper bound
}

transformed data {
  vector[N] ycen;                             // centered response
  real meany;
  vector[K] mustartcen;
  
  meany <- mean(y);
  ycen <- y-meany;
  mustartcen <- mustart-meany;
  
}

parameters {
  simplex[K] theta;                           // mixing proportions; simplex will sum to 1
  ordered[K] mu;                              // locations of mixture components
  real<lower=0, upper=sigmaupper> sigma[K];   // scales of mixture components
}

transformed parameters{

}

model {
  // model calculations
  real ps[K];                                 // initial log component densities
  
  for (n in 1:N){
    for (k in 1:K){
      ps[k] <- log(theta[k]) + normal_log(ycen[n], mu[k], sigma[k]);
    }
  }

  // priors
  theta ~ beta(1,1);
  //sigma ~ cauchy(0, 2.5);
  sigma ~ cauchy(0, 5);

  for (k in 1:K){
    mu[k] ~ normal(mustartcen[k], sigmastart);
    //mu[k] ~ normal(0, 1);
  }

  // likelihood
  increment_log_prob(log_sum_exp(ps));
}

generated quantities {
  vector[K] finalmu;
  for (k in 1:K){
    finalmu[k] <-  mu[k]+meany;                // back to original scale
  }
}
'



### Test Run ###
library(rstan)
test <- stan(model_code = stanmodelcode, model_name = "example", 
            data = standat, iter = 7000, warmup=2000, thin=5, chains = 2, 
            verbose = F)   

traceplot(test)
print(test, digits=3)


library(flexmix)
flexmod1 = flexmix(eruptions~1, k=2, data=faithful, control=list(tolerance=1e-12, iter.max=1000))
summary(flexmod1)
parameters(flexmod1)

flexmod2 = flexmix(waiting~1, k=2, data=faithful, control=list(tolerance=1e-8, iter.max=1000))
summary(flexmod2)
parameters(flexmod2)

flexmod3 = flexmix(y~1, k=2, control=list(tolerance=1e-8, iter.max=1000))
summary(flexmod3)
parameters(flexmod3)

### try hard
library(parallel)
cl = makeCluster(3)
clusterEvalQ(cl, library(rstan))

clusterExport(cl, c('stanmodelcode', 'standat', 'test')) 

p = proc.time()
parfit = parSapply(cl, 1:3, function(i) stan(model_code = stanmodelcode, model_name = "mixture", fit = test, 
                                             data = standat, iter = 62000, warmup=12000, thin=50, chains = 1, chain_id=i,
                                             verbose = T), 
                   simplify=F) 

proc.time() - p


fit = sflist2stanfit(parfit) 
print(fit, digits=3)
traceplot(fit)
# pairs(fit)

mus = extract(fit, par='finalmu')$finalmu
library(ggplot2); library(reshape2)
gdat = melt(mus)

ggplot(aes(x=waiting), data=faithful) +
  geom_vline(aes(xintercept=value), col='gray90',alpha=.05, data=gdat) +
  geom_density() +
  geom_point(aes(x=value, y=0), alpha=.05, data=gdat) +
  ggtheme

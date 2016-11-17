#---------------------------------------------------------------------------------#
# Mixture model (draft); The Stan manual notes issues with mixture modeling, esp. #
# label switching (as does the bugs book), but the code in the manual chapter on  #
# mixture modeling does not appear able to reproduce the old faithful mixture or  #
# even a simulated mixture of normals, even with some heavy hand holding on       # 
# starting points, limits etc.  On occasion the means can be recovered, but the   #
# variances are all over the place and notably influenced by the prior.  The      #
# Stan manual chapter on problematic posteriors covers many of the issues.        #
#                                                                                 #
# The following stan code however is verbatim from github and should work fine,   #
# albeit slowly.                                                                  #
# github stan/src/models/basic_estimators/normal_mixture_k_prop.stan              #
#---------------------------------------------------------------------------------#

data(faithful)
head(faithful)
par(mfrow=c(1,2))
apply(faithful, 2, ggplot2:::qplot, geom='density')
layout(1)


y = rnorm(500, c(50,80), c(5,10))
ggplot2:::qplot(y, geom='density')
psych::describe(data.frame(y[seq(1,500, 2)], y[seq(2,500, 2)]))


# take your pick
standat = list(N=nrow(faithful), K=2, y=faithful$waiting)

standat = list(N=nrow(faithful), K=2, y=faithful$eruptions)

standat = list(N=length(y), K=2, y=y)


stanmodelcode = '
data {
  int<lower=1> K;                                                  # K components
  int<lower=1> N;                                                  # N observations
  real y[N];                                                       # variable of interest
}

parameters {
  simplex[K] theta;                                                # mixing proportions
  simplex[K] mu_prop;
  real mu_loc;
  real<lower=0> mu_scale;
  real<lower=0> sigma[K];                                          # sds of the components
}

transformed parameters {
  ordered[K] mu;
  mu = mu_loc + mu_scale * cumulative_sum(mu_prop);               # means of the components
}

model {
  // prior
  mu_loc ~ cauchy(0,5);
  mu_scale ~ cauchy(0,5);
  sigma ~ cauchy(0,5);

  // likelihood
  {
    real ps[K];
    vector[K] log_theta;
    log_theta = log(theta);

    for (n in 1:N) {
      for (k in 1:K) {
        ps[k] = log_theta[k] + normal_lpdf(y[n] | mu[k], sigma[k]);
      }
      target += log_sum_exp(ps);
    }
  }
}
'


################
### Test Run ###
################
library(rstan)

# the following may take several minutes per chain depending on the data
test = stan(model_code = stanmodelcode, model_name = "example", 
            data = standat, iter = 7000, warmup=2000, thin=5, chains = 2, cores=2,
            verbose = F)

# shinystan::launch_shinystan(test)
print(test, digits=3)

### compare to flexmix: coefs in flexmod1 = mu in test
library(flexmix)
flexmod1 = flexmix(standat$y~1, k=2, control=list(tolerance=1e-12, iter.max=1000))
summary(flexmod1)
parameters(flexmod1)


######################
### Production run ###
######################

# This can take a notably long time (i.e. hours) depending on the data.
fit = stan(model_code = stanmodelcode, model_name = "mixture", fit = test, 
                                             data = standat, iter = 62000, warmup=12000, thin=50, cores = 4,
                                             verbose = T)
print(fit, digits=3)
shinystan::launch_shinystan(fit)
 

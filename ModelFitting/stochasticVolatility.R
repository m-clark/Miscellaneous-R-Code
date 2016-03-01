### Stochastic Volatility Model for centered time series over t equally spaced
### points. latent parameter h is the log volatility, phi the persistence of the
### volatility and mu the mean log volatility. ϵ is the white-noise shock and δ
### the shock on volatility.  The stan code is based on that in the manual.
# y_t = exp(h_t/2)*ϵ_t
# h_t = μ + φ*(h_{t-1}-μ) + δ_t*σ
# h_1 ~ N(μ, σ/sqrt(1-φ^2))
# ϵ_t ~ N(0,1); δ_t ~ N(0,1) 
# 
# With rearranging,
# ϵ_t = y_t*exp(-h_t/2)
# y_t ~ N(0, exp(h_t/2)
# h_t ~ N(μ + φ*(h_t-μ), σ)

### data available in data repository; regards inflation based on the U.S. consumer 
### price index (infl =  400*log(cpi_t/cpi_{t-1})), from the second quarter of 
### 1947 to the second quarter of 2011 (from Statistical Computation and
### Modeling 2014, chap 11)

d = read.csv('../Datasets/us cpi/USCPI.csv', header=F)
inflation = d[,1]
summary(inflation)
inflation_cen = scale(inflation, scale=F)


# this original code keeps to the above formulation but can take a long time to converge. 
# ϵ_t and δ_t are implicit
stanmodelcodeConceptual = "
data {
  int<lower=0> T;                      // # time points (equally spaced)
  vector[T] y;                         // mean corrected response at time t
}

parameters {
  real mu;                             // mean log volatility
  real<lower=-1,upper=1> phi;          // persistence of volatility
  real<lower=0> sigma;                 // white noise shock scale
  vector[T] h;                         // log volatility at time t
}

model {
  //priors
  phi ~ uniform(-1,1);
  sigma ~ cauchy(0,5);
  mu ~ cauchy(0,10);

  //likelihood
  h[1] ~ normal(mu, sigma / sqrt(1 - phi * phi));
  for (t in 2:T)
    h[t] ~ normal(mu + phi * (h[t - 1] - mu), sigma);

  for (t in 1:T)
    y ~ normal(0, exp(h[t] / 2));
}
"


stanmodelcodePerformance = "
data {
  int<lower=0> T;                      // time points (equally spaced)
  vector[T] y;                         // mean corrected response at time t
}

parameters {
  real mu;                             // mean log volatility
  real<lower=-1,upper=1> phi;          // persistence of volatility
  real<lower=0> sigma;                 // white noise shock scale
  vector[T] h_std;                     // standardized log volatility at time t
}

transformed parameters{
  vector[T] h;                         // log volatility at time t
  h <- h_std * sigma;
  h[1] <- h[1] / sqrt(1-phi * phi);
  h <- h + mu;
  for (t in 2:T)
    h[t] <- h[t] + phi * (h[t-1] - mu);
}

model {
  //priors
  phi ~ uniform(-1,1);
  sigma ~ cauchy(0,5);
  mu ~ cauchy(0,10);
  h_std ~ normal(0,1);

  //likelihood
  y ~ normal(0, exp(h/2));
}

generated quantities{
  vector[T] yRep;

  for (t in 1:T){
    yRep[t] <- normal_rng(0, exp(h[t]/2));
  }  
  
}
"

# use c() to get rid of matrix format or specify as matrix instead of vector in model code
standat = list(T=length(inflation_cen), y=c(inflation_cen))

library(rstan)

### test runs
# # very slow
# stochvol_1 = stan(model_code = stanmodelcodeConceptual,
#                  data = standat, chains = 2,
#                  iter = 4000, warmup = 2000, thin = 10)

# print(stochvol_1, dig=3, par=c('mu', 'phi', 'sigma'), probs=c(.025,.5,.975))
# traceplot(stochvol_1,  par=c('mu', 'phi', 'sigma'))

# # much better
stochvol_2 = stan(model_code = stanmodelcodePerformance,
                 data = standat, chains = 2,
                 iter = 4000, warmup = 2000, thin = 10)

# print(stochvol_2, dig=3, par=c('mu', 'phi', 'sigma'), probs=c(.025,.5,.975))
# traceplot(stochvol_2,  par=c('mu', 'phi', 'sigma'))


library(parallel)
cl = makeCluster(4)
clusterEvalQ(cl, library(rstan))
clusterExport(cl, c('standat', 'stanmodelcodePerformance'))

stochvol = parSapply(cl, 1:4, function(chain) stan(model_code = stanmodelcodePerformance,
                data = standat, chains = chain,
                iter = 22000, warmup = 2000, thin = 20))
stopCluster(cl)

stochvol_final = sflist2stanfit(stochvol)
print(stochvol_final, dig=3, par=c('mu', 'phi', 'sigma'), probs=c(.025,.5,.975))

traceplot(stochvol_final,  par=c('mu', 'phi', 'sigma'))
# pairs(stochvol_final,  par=c('mu', 'phi', 'sigma'))
plot(stochvol_final)


# get log volatility and data replicates
h = extract(stochvol_final, 'h')$h
yRep = apply(h, 1, function(h) rnorm(length(inflation), 0, exp(h/2)))
yRep2 = extract(stochvol_final, 'yRep')$yRep
dim(yRep)


library(lubridate); library(scales)
series = ymd(paste0(rep(1947:2014, e=4),'-', c('01','04','07','10'), '-', '01'))
seriestext = series[1:length(inflation)]

# compare to fig. 11.1 in the text
plot(seriestext, colMeans(h), cex=.5, col='gray50', ylim=c(-2,4), type='l', bty='n')

plot(seriestext, inflation_cen, cex=.5, col='gray50', ylim=c(-10,20), type='l', bty='n')
sapply(sample(1:dim(yRep)[2], 100), function(i) lines(seriestext, yRep[,i], cex=.5, col=alpha('#FF5500', .05)) )
lines(seriestext, inflation_cen, cex=.5, col='gray50', ylim=c(-10,20), type='l', bty='n')


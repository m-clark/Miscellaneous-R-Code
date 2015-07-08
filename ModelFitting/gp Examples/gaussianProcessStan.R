# Data and parameter setup ------------------------------------------------

# Data
set.seed(1234)
N = 20
Ntest = 200
x = rnorm(N, sd=5)
y = sin(x) + rnorm(N, sd=.1)
xtest = seq(min(x)-1, max(x)+1, l=Ntest)
plot(x,y, pch=19, col='#ff5500')

# parameters
eta_sq = 1
rho_sq = 1
sigma_sq = .1

# Covariance function same as implemented in the Stan code.
Kfn <- function (x, eta_sq, rho_sq, sigma_sq) {
  N = length(x)
  Sigma = matrix(NA, N, N)
  
  # off diag elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * (x[i] - x[j])^2);
      Sigma[j,i] <- Sigma[i,j];
    }
  }
  
  # diagonal elements
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq; # + jitter
  Sigma
}



# Vis prior functions -----------------------------------------------------
xinit = seq(-5,5,.2)
xprior = MASS::mvrnorm(3, 
                       mu=rep(0, length(xinit)), 
                       Sigma=Kfn(x=xinit,
                                 eta_sq = eta_sq, 
                                 rho_sq = rho_sq, 
                                 sigma_sq = sigma_sq))


library(reshape2)
gdat = melt(data.frame(x=xinit, y=t(xprior)), id='x')

library(ggvis)
gdat %>% 
  ggvis(~x, ~value) %>% 
  group_by(variable) %>% 
  layer_paths(strokeOpacity:=.5) %>% 
  add_axis('x', grid=F) %>% 
  add_axis('y', grid=F)



# Stan model code ---------------------------------------------------------

gp = "
data {
  int<lower=1> N;                                # initial sample size
  vector[N] x;                                   # covariate
  vector[N] y;                                   # target
  int<lower=0> Ntest;                            # prediction set sample size
  vector[Ntest] xtest;                           # prediction values for covariate
}

transformed data {
  vector[N] mu;
  
  mu <- rep_vector(0, N);                        # mean function
}

parameters {
  real<lower=0> eta_sq;                          # parameters of squared exponential covariance function
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
}

transformed parameters {
  real<lower=0> rho_sq;
  rho_sq <- inv(inv_rho_sq);
}

model {
  matrix[N,N] Sigma;

  # off-diagonal elements for covariance matrix
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * pow(x[i] - x[j],2));
      Sigma[j,i] <- Sigma[i,j];
    }
  }

  # diagonal elements
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq;             # + jitter for pos def

  # priors
  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);

  # sampling distribution
  y ~ multi_normal(mu,Sigma);
}

generated quantities {
  vector[Ntest] muTest;                          # The following produces the posterior predictive draws
  vector[Ntest] yRep;                            # see GP section of Stan man- 'Analytical Form...'
  matrix[Ntest,Ntest] L;
  {
  matrix[N,N] Sigma;
  matrix[Ntest,Ntest] Omega;
  matrix[N,Ntest] K;
  matrix[Ntest,N] K_transpose_div_Sigma;
  matrix[Ntest,Ntest] Tau;
  
  # Sigma
  for (i in 1:N)
  for (j in 1:N)
  Sigma[i,j] <- eta_sq * exp(-pow(x[i] - x[j], 2)) + if_else(i==j, sigma_sq, 0.0);
  
  # Omega
  for (i in 1:Ntest)
  for (j in 1:Ntest)
  Omega[i,j] <- eta_sq * exp(-pow(xtest[i] - xtest[j], 2)) + if_else(i==j, sigma_sq, 0.0);
  
  # K
  for (i in 1:N)
  for (j in 1:Ntest)
  K[i,j] <- eta_sq * exp(-pow(x[i] - xtest[j], 2));
  
  K_transpose_div_Sigma <- K' / Sigma;
    muTest <- K_transpose_div_Sigma * y;
    Tau <- Omega - K_transpose_div_Sigma * K;
    
    for (i in 1:(Ntest-1))
    for (j in (i+1):Ntest)
    Tau[i,j] <- Tau[j,i];
    
    L <- cholesky_decompose(Tau);
  }
          
  yRep <- multi_normal_cholesky_rng(muTest, L);
}
"



# Compile Check -----------------------------------------------------------

standata = list(N=N, x=x, y=y, xtest=xtest, Ntest=200)

library(rstan)
fit0 = stan(data=standata, model_code = gp, iter = 1, chains=1)



# Main Run ----------------------------------------------------------------

iterations = 12000
wu = 2000
th = 20
chains = 4


# fit = stan(data=standata, model_code = gp, iter = 12000, warmup = 2000, thin=20, 
#             chains=4, fit = fit0)
# 
# fit


# With N = 20 Ntest = 200 takes about 2 min
library(parallel)
cl = makeCluster(4)
clusterExport(cl, c('gp',  'standata', 'fit0','iterations', 'wu', 'th', 'chains'))
clusterEvalQ(cl, library(rstan))


p = proc.time()
fit = parLapply(cl, seq(chains), function(chain) stan(data=standata, model_code = gp, 
                                                      iter = iterations,
                                                      warmup = wu, thin=th, 
                                                      chains=1, chain_id=chain, 
                                                      fit = fit0)
                )
(proc.time() - p)/3600

stopCluster(cl)



# Summarize and Vis -------------------------------------------------------

fit = sflist2stanfit(fit)
# takes a bit to print
# print(fit, par=c('eta_sq','rho_sq','sigma_sq'))

# library(shinyStan)
# launch_shinystan(fit)


# Extract and visualize posterior predictive draws
yRep = extract(fit, 'yRep')$yRep

gdat = data.frame(x,y)
gdat2 = melt(data.frame(x = sort(xtest), y=t(yRep[sample(2000, 3),])), id='x')

gdat2 %>% 
  ggvis(~x, ~value) %>% 
  group_by(variable) %>% 
  layer_paths(strokeOpacity:=.25) %>% 
  layer_points(x=~x, y=~y, fill:='#ff5500', data=gdat) %>% 
  add_axis('x', grid=F) %>% 
  add_axis('y', grid=F)  


# Visualize fit
yRepMean = get_posterior_mean(fit, 'yRep')[,5]
gdat3 = data.frame(x = sort(xtest), y=yRepMean)

gdat3 %>% 
  ggvis(~x, ~y) %>% 
  layer_paths(strokeOpacity:=.5, stroke:='blue') %>% 
  layer_points(x=~x, y=~y, fill:='#ff5500', data=gdat) %>% 
  add_axis('x', grid=F) %>% 
  add_axis('y', grid=F)  



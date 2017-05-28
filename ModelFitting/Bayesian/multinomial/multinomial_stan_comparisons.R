stan_catlogit_km1 = "
data {
  int<lower=1> K;                   # number of classes
  int<lower=1> N;                   # nrow of x
  int<lower=1> D;                   # ncol of x
  int<lower=1, upper=K> y[N];       # target as integer
  vector[D] x[N];                   # array of D
}

transformed data {
  row_vector[D] zeros;              # create reference level coefs of zero
  zeros = rep_row_vector(0, D);
}

parameters {
  matrix[K-1,D] beta_raw;           # estimated coefs
}

transformed parameters{
  matrix[K, D] beta;
  beta = append_row(zeros, beta_raw);
}

model {
  # prior
  to_vector(beta_raw) ~ normal(0, 10);
  
  # likelihood
  for (n in 1:N)
    y[n] ~ categorical_logit(beta * x[n]);
}"


# Applies some vectorization, some speed gain; but looks negligible
stan_catlogit_km1_vec = "
data {
  int<lower=1> K;                   # number of classes
  int<lower=1> N;                   # nrow of x
  int<lower=1> D;                   # ncol of x
  int<lower=1, upper=K> y[N];       # target as integer
  matrix[N,D] x;                    # model matrix
}

transformed data{
  matrix[D,N] xt;
  row_vector[D] zeros;
  
  xt = x';
  zeros = rep_row_vector(0, D);
}

parameters {
  matrix[K-1,D] beta_raw;
}

transformed parameters{
  matrix[K, D] beta;
  beta = append_row(zeros, beta_raw);
}

model {
  matrix[K,N] L;                    # Linear predictor

  L = beta * xt;


  to_vector(beta_raw) ~ normal(0, 10);
  
  for (n in 1:N)
    y[n] ~ categorical_logit(L[,n]);
}
"

# keeps things in columnar format, no data transpose
stan_catlogit_km1_vec_mlogit = "
data {
  int K;
  int N;
  int D;
  int y[N];
  matrix[N,D] x;
}

transformed data{
  vector[D] zeros;
  
  zeros = rep_vector(0, D);
}

parameters {
  matrix[D,K-1] beta_raw;
}

transformed parameters{
  matrix[D, K] beta;
  beta = append_col(zeros, beta_raw);
}

model {
  matrix[N, K] L;                   # Linear predictor

  L = x * beta;

  to_vector(beta_raw) ~ normal(0, 10);
  
  for (n in 1:N)
    y[n] ~ categorical_logit(to_vector(L[n]));
}
"


# somewhat slower
stan_catlogit_km1_vec_constraint = "
data {
  int K;
  int N;
  int D;
  int y[N];
  matrix[N,D] x;
}

transformed data{
  matrix[D,N] xt;

  xt = x';
}

parameters {
  matrix[K-1,D] beta_raw;

}

transformed parameters{
  matrix[K, D] beta;

  for (d in 1:D) {
    beta[1,d] = -sum(beta_raw[,d]);
    beta[2:K,d] = beta_raw[,d];
  }
}

model {
  matrix[K,N] L;

  L = beta * xt;                    # Linear predictor
  
  to_vector(beta_raw) ~ normal(0, 10);

  for (n in 1:N)
    y[n] ~ categorical_logit(L[,n]);
}
"
 
# this follows the stan manual example, but is slow, adds to code verbosity, and
# doesn't provide correct (though close) estimates. The scale parameters in
# particular are problematic.
# stan_catlogit_km1_vec_simplex = "
# data {
#   int K;
#   int N;
#   int D;
#   int y[N];
#   matrix[N,D] x;
# }
# 
# transformed data{
#   matrix[D,N] xt;
#   row_vector[D] zeros;
#   
#   xt = x';
# }
# 
# parameters {
#   simplex[K] beta_raw1;
#   simplex[K] beta_raw2;
#   simplex[K] beta_raw3;
#   simplex[K] beta_raw4;
# 
#   vector<lower=0>[4] beta_scale;
# }
# 
# transformed parameters{
#   matrix[K, D] beta;
# 
# 
#   beta[,1] = beta_scale[1] * (beta_raw1 - 1.0 / K);
#   beta[,2] = beta_scale[2] * (beta_raw2 - 1.0 / K);
#   beta[,3] = beta_scale[3] * (beta_raw3 - 1.0 / K);
#   beta[,4] = beta_scale[4] * (beta_raw4 - 1.0 / K);
# 
# }
# 
# model {
#   matrix[K,N] L;
# 
#   L = beta * xt;
#   
#   beta_raw1 ~ normal(0, 10);
#   beta_raw2 ~ normal(0, 10);
#   beta_raw3 ~ normal(0, 10);
#   beta_raw4 ~ normal(0, 10);
# 
#   # beta_scale ~ student_t(3, 0, 100);
# 
#   for (n in 1:N)
#     y[n] ~ categorical_logit(L[,n]);
# }
# "





# Import data and setup  --------------------------------------------------

library(haven); library(tidyverse)
program = read_dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta") %>% 
  as_factor() %>% 
  mutate(prog = relevel(prog, ref = "academic"))

library(mlogit)
head(program[,1:5])
programLong = program %>% 
  select(id, prog, ses, write) %>% 
  mlogit.data(data=, shape='wide', choice='prog', id.var='id')
head(programLong)
mlogit_mod = mlogit(prog ~ 1|ses + write, data=programLong)
mlogit_coefs = coef(mlogit_mod)[c(1,3,5,7,2,4,6,8)]

X = model.matrix(prog ~ ses + write, data = program)
y = program$prog
X = X[order(y),]
y = y[order(y)]

# N = sample size, x is the model matrix, y integer version of class outcome, k=
# number of classes, D is dimension of model matrix
datalist = list(N=nrow(X), x = X, y=as.integer(y), K=n_distinct(y), D=ncol(X))


library(rstan)

# categorical based on K-1
bayes_catlogit_km1 = stan(model_code=stan_catlogit_km1, data=datalist, cores=4)
bayespar_catlogit_km1 = get_posterior_mean(bayes_catlogit_km1, par='beta_raw')[,5]
cbind(mlogit_coefs, bayespar_catlogit_km1)


# categorical based on K-1 + vec
bayes_catlogit_km1_vec = stan(model_code=stan_catlogit_km1_vec, data=datalist, cores=4)
bayespar_catlogit_km1_vec = get_posterior_mean(bayes_catlogit_km1_vec, par='beta_raw')[,5]
cbind(mlogit_coefs, bayespar_catlogit_km1_vec)

bayes_catlogit_km1_vec_mlogit = stan(model_code=stan_catlogit_km1_vec_mlogit, data=datalist, cores=4)
bayespar_catlogit_km1_vec_mlogit = get_posterior_mean(bayes_catlogit_km1_vec_mlogit, par='beta_raw')[,5]
cbind(coef(mlogit_mod), bayespar_catlogit_km1_vec_mlogit)


# categorical based on K-1 + vec constraint
bayes_catlogit_km1_vec_constraint = stan(model_code=stan_catlogit_km1_vec_constraint, data=datalist, cores=4)
bayespar_catlogit_km1_vec_constraint = get_posterior_mean(bayes_catlogit_km1_vec_constraint, par='beta')[,5]
bayespar_catlogit_km1_vec_constraint = c(bayespar_catlogit_km1_vec_constraint[5:12]-bayespar_catlogit_km1_vec_constraint[1:4])  
cbind(mlogit_coefs, bayespar_catlogit_km1_vec_constraint)




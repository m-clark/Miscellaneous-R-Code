data {                                           
  int<lower=1> N;                                # sample size
  int <lower=1> Ncol_x;                          # number of columns of teh covariate matrix
  int<lower=1> I;                                # number of groups
  vector[N] y;                                   # target variable: center or standardize
  matrix[N, Ncol_x] X;                           # covariate matrix: center or standardize; currently set up for one covariate
  int<lower=1,upper=I> Group[N];                 # grouping factor
}

transformed data {
  int <lower=1> P;
//   matrix[N,Ncol_x+1] MM;                      # in case you want to work with a model matrix
//   
//   MM = append_col(rep_vector(1, N), X);
//   P = cols(MM);
  P = Ncol_x + 1;
}

parameters {
  vector[P] fixefs;                              # fixed effects
  real<lower=0> sigma_int_fe;                    # fixed effects scales
  real<lower=0> sigma_beta_fe;                      
  cholesky_factor_corr[P] Omega_FE_chol;         # correlation matrix for fixed effects (chol decomp)
  
  real<lower=0> sigma_int_re;                    # random effects scales
  real<lower=0> sigma_beta_re;
  real<lower=0> sigma_int_re_hyper;              # random effects scales hyperprior
  real<lower=0> sigma_beta_re_hyper;

  vector[2] gamma[I];                            # individual effects; array allows vectorized multinormal
  cholesky_factor_corr[2] Omega_RE_chol;         # correlation matrix for random intercepts and slopes (chol decomp)
  
  real<lower=0> sigma_y;                         # residual sd
}

transformed parameters {
  real Intercept;
  real beta;

  Intercept = fixefs[1];
  beta = fixefs[2];
} 

model {
  matrix[P,P] DC_fe;
  matrix[2,2] DC_re;
  vector[P] sigma_fe;
  vector[2] sigma_re;
  vector[N] LP;                                  # Linear predictor

  # priors
  Omega_FE_chol ~  lkj_corr_cholesky(2.0); 
  Omega_RE_chol ~  lkj_corr_cholesky(2.0); 
  
  ### hyperpriors
  sigma_int_re_hyper ~ exponential(.2);
  sigma_beta_re_hyper ~ exponential(.2);

  ### priors
  sigma_int_fe ~ cauchy(0, 2.5);        
  sigma_beta_fe ~ cauchy(0, 2.5); 
  
  sigma_int_re ~ cauchy(0, sigma_int_re_hyper);        
  sigma_beta_re ~ cauchy(0, sigma_beta_re_hyper); 

  sigma_y ~ cauchy(0, 2.5);                      
  
  # multivariate draw for fixed effects
  sigma_fe[1] = sigma_int_fe;
  sigma_fe[2] = sigma_beta_fe;
  DC_fe = diag_pre_multiply(sigma_fe, Omega_FE_chol);

  fixefs ~ multi_normal_cholesky(rep_vector(0, 2), DC_fe);
  
  
  # multivariate draw for random effects
  sigma_re[1] = sigma_int_re;
  sigma_re[2] = sigma_beta_re;
  DC_re = diag_pre_multiply(sigma_re, Omega_RE_chol);
  
  gamma ~ multi_normal_cholesky(fixefs, DC_re);
  
  
  # likelihood
  for (n in 1:N){                          
    # LP[n] = MM[n] * betavec;                  # tried matrix approach, but it was actually 2-3 times slower
    LP[n] = gamma[Group[n],1] + gamma[Group[n],2] * X[n,1];
  }
  y ~ normal(LP, sigma_y);
}

generated quantities {
  corr_matrix[P] Omega_FE;                       # correlation of FE
  corr_matrix[2] Omega_RE;                       # correlation of RE

  Omega_FE = tcrossprod(Omega_FE_chol);
  Omega_RE = tcrossprod(Omega_RE_chol);
}
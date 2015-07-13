data {
  int<lower=1> N;                                # initial sample size
  int<lower=1> D;                                # number of covariates
  matrix[N,D] X;                                 # covariate matrix
  vector[N] y;                                   # target
  int<lower=0> Ntest;                            # prediction set sample size
  matrix[Ntest,D] Xtest;                           # prediction values for covariate
  int<lower=1> K;                                # number of factors
}

transformed data {
  vector[N] mu;

  mu <- rep_vector(0, N);                        # mean function
}

parameters {
  real<lower=0, upper=1> eta_sq;                 # parameters of squared exponential covariance function
  real<lower=0, upper=1> sigma_sq;               # eta_sq + sigma_sq = var of predicted y; eta_sq is variance explained by the function so can put an upper limit
  vector<lower=0, upper=100>[D] l_sq;               # characteristic length
  real<lower=-2, upper=2> lambda;                # factor loadings, for this problem a single value, as one of them will be fixed to 1 for identification purposes

//   real<lower=0> eta_sq_scale;                    # hyperpriors
//   real<lower=0> sigma_sq_scale;
//   real<lower=0> l_sq_scale;
}

transformed parameters {
  vector[D] inv_l_sq;
  matrix[D,K] L;

  L[1,1] <- 1;
  L[2,1] <- lambda;

  for (d in 1:D) inv_l_sq[d] <- 1/l_sq[d];
}

model {
  matrix[N,N] Sigma;


  # off-diagonal elements for covariance matrix
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * exp(-.5 * quad_form(tcrossprod(L) + diag_matrix(inv_l_sq), to_vector(X[i] - X[j])));
      Sigma[j,i] <- Sigma[i,j];
    }
  }

  # diagonal elements
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq;             # + jitter for pos def

  # hyperpriors
//   eta_sq_scale ~ exponential(.2);
//   sigma_sq_scale ~ exponential(.2);
//   l_sq_scale ~ exponential(.2);

  # priors
  eta_sq ~ cauchy(0, 5);
  sigma_sq ~ cauchy(0, 5);
  l_sq ~ cauchy(0, 5);
  lambda ~ normal(0, 1);

  # sampling distribution
  y ~ multi_normal(mu, Sigma);
}

generated quantities {
  vector[Ntest] muTest;                          # The following produces the posterior predictive draws
  vector[Ntest] yRep;                            # see GP section of Stan man- 'Analytical Form...'
  matrix[Ntest,Ntest] Q;
  {
    matrix[N,N] Sigma;
    matrix[Ntest,Ntest] Omega;
    matrix[N,Ntest] K_;
    matrix[Ntest,N] K_transpose_div_Sigma;
    matrix[Ntest,Ntest] Tau;

    # Sigma
    for (i in 1:N)
      for (j in 1:N)
        Sigma[i,j] <- eta_sq * exp(-.5 * quad_form(tcrossprod(L) + diag_matrix(inv_l_sq), to_vector(X[i] - X[j]))) + if_else(i==j, sigma_sq, 0.0);

    # Omega
    for (i in 1:Ntest)
      for (j in 1:Ntest)
        Omega[i,j] <- eta_sq * exp(-.5 * quad_form(tcrossprod(L) + diag_matrix(inv_l_sq), to_vector(Xtest[i] - Xtest[j]))) + if_else(i==j, sigma_sq, 0.0);

    # K
    for (i in 1:N)
      for (j in 1:Ntest)
        K_[i,j] <- eta_sq * exp(-.5 * quad_form(tcrossprod(L) + diag_matrix(inv_l_sq), to_vector(X[i] - Xtest[j]))); 

    K_transpose_div_Sigma <- K_' / Sigma;
    muTest <- K_transpose_div_Sigma * y;
    Tau <- Omega - K_transpose_div_Sigma * K_;

    for (i in 1:(Ntest-1))
      for (j in (i+1):Ntest)
        Tau[i,j] <- Tau[j,i];

    Q <- cholesky_decompose(Tau);
  }

  yRep <- multi_normal_cholesky_rng(muTest, Q);
}
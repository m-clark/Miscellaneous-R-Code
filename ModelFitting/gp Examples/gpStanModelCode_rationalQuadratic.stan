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
  real<lower=0> alpha;                           # parameters of rational quadratic covariance function
  real<lower=0> l_;
  real<lower=0> eta_sq;
  real<lower=0> sigma_sq;

  real<lower=0> alpha_scale;                     # hyperpriors
  real<lower=0> l_scale;
  real<lower=0> eta_sq_scale;
  real<lower=0> sigma_sq_scale;
}

transformed parameters {
}

model {
  matrix[N,N] Sigma;

  # off-diagonal elements for covariance matrix
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * (1 + (1/(2*alpha)) * pow(fabs(x[i] - x[j]) /l_, 2))^-alpha;
      Sigma[j,i] <- Sigma[i,j];
    }
  }

  # diagonal elements
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq;             # + jitter for pos def

  # hyperpriors
  eta_sq_scale ~ exponential(.2);
  l_scale ~ exponential(.2);
  sigma_sq_scale ~ exponential(.2);
  alpha_scale ~ exponential(.2);

  # priors
  eta_sq ~ cauchy(0, eta_sq_scale);
  l_ ~ cauchy(0, l_scale);
  sigma_sq ~ cauchy(0, sigma_sq_scale);
  alpha ~ cauchy(0, alpha_scale);

  # sampling distribution
  y ~ multi_normal(mu, Sigma);
}

generated quantities {
  vector[Ntest] muTest;                          # The following produces the posterior predictive draws
  vector[Ntest] yRep;                            # see GP section of Stan man- 'Analytical Form...'
  matrix[Ntest,Ntest] L;
  {
    matrix[N, N] Sigma;
    matrix[Ntest, Ntest] Omega;
    matrix[N,Ntest] K;
    matrix[Ntest,N] K_transpose_div_Sigma;
    matrix[Ntest, Ntest] Tau;

    # Sigma
    for (i in 1:N)
      for (j in 1:N)
        Sigma[i,j] <- eta_sq * (1 + (1/(2*alpha)) * pow(fabs(x[i] - x[j]) /l_, 2))^-alpha + if_else(i==j, sigma_sq, 0.0);

    # Omega
    for (i in 1:Ntest)
      for (j in 1:Ntest)
        Omega[i,j] <- eta_sq * (1 + (1/(2*alpha)) * pow(fabs(xtest[i] - xtest[j]) /l_, 2))^-alpha + if_else(i==j, sigma_sq, 0.0);

    # K
    for (i in 1:N)
      for (j in 1:Ntest)
        K[i,j] <- eta_sq * (1 + (1/(2*alpha)) * pow(fabs(x[i] - xtest[j]) /l_, 2))^-alpha;

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
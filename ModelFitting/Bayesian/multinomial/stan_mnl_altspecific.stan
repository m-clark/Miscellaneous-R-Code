data {
  int K;                               # number of choices
  int N;                               # number of individuals
  int D;                               # number of indiv specific vars
  int G;                               # number of alt specific vars
  int T;                               # number of alt constant vars
  
  int y[N*K];                          # choices
  vector[N*K] choice;                  # choice made (logical)
  
  matrix[N,D] X;                       # data for indiv specific effects
  matrix[N*K, G] Y;                    # data for alt specific effects
  matrix[N*(K-1),T] Z;                 # data for alt constant effects
  
}

transformed data {

}

parameters {   
  matrix[D, K-1] beta;                 # individual specific coefs
  matrix[G, K]  gamma;                 # choice specific coefs for alt-specific variables
  vector[T]  theta;                    # choice constant coefs for alt-specific variables
  
}

transformed parameters{
  
}

model {
  matrix[N, K-1] Vx;                   # Utility for individual vars
  
  vector[N*K] Vy0;
  matrix[N, K-1] Vy;                   # Utility for alt-specific/alt-varying vars
  
  vector[N*(K-1)] Vz0;
  matrix[N, (K-1)] Vz;                 # Utility for alt-specific/alt-constant vars

  matrix[N,K-1] V;                     # combined utilities
  
  vector[N] baseProbVec;               # reference group probabilities
  real ll0;                            # intermediate log likelihood
  real loglik;                         # final log likelihood

  
  to_vector(beta) ~ normal(0, 10);     # diffuse priors on coefficients
  to_vector(gamma) ~ normal(0, 10);    # diffuse priors on coefficients
  to_vector(theta) ~ normal(0, 10);    # diffuse priors on coefficients
  

  Vx = X * beta;
  
  for(alt in 1:K){
    vector[G] par;
    int start;
    int end;

    par = gamma[,alt];
    start = N*alt-N+1;
    end = N*alt;
    Vy0[start:end] = Y[start:end,] * par;
    if(alt>1) Vy[,alt-1] = Vy0[start:end] - Vy0[1:N];
  }
  
  Vz0 = Z * theta;
  
  for(alt in 1:(K-1)){
    int start;
    int end;

    start = N*alt-N+1;
    end = N*alt;
    Vz[,alt] = Vz0[start:end];
  }

  V = Vx + Vy + Vz;

  for(n in 1:N)  baseProbVec[n] = 1/(1 + sum(exp(V[n])));
  ll0 = dot_product(to_vector(V), choice[(N+1):(N*K)]); # just going to assume no neg index
  loglik = sum(log(baseProbVec)) + ll0;
  target += loglik;
  
}


generated quantities {
  matrix[N,K-1] fitted_nonref;
  vector[N] fitted_ref;
  matrix[N,K] fitted;
  
  matrix[N,K-1] Vx;                    # Utility for individual vars
  
  vector[N*K] Vy0;
  matrix[N,K-1] Vy;                    # Utility for alt-specific/alt-varying vars
  
  vector[N*(K-1)] Vz0;
  matrix[N, (K-1)] Vz;                 # Utility for alt-specific/alt-constant vars

  matrix[N,K-1] V;                     # combined utilities
  
  vector[N] baseProbVec;               # reference group probabilities

  Vx = X * beta;
  
  for(alt in 1:K){
    vector[G] par;
    int start;
    int end;

    par = gamma[,alt];
    start = N*alt-N+1;
    end = N*alt;
    Vy0[start:end] = Y[start:end,] * par;
    if(alt>1) Vy[,alt-1] = Vy0[start:end] - Vy0[1:N];
  }
  
  Vz0 = Z * theta;
  
  for(alt in 1:(K-1)){
    int start;
    int end;

    start = N*alt-N+1;
    end = N*alt;
    Vz[,alt] = Vz0[start:end];
  }

  V = Vx + Vy + Vz;
  
  for(n in 1:N)  baseProbVec[n] = 1/(1 + sum(exp(V[n])));
  fitted_nonref = exp(V) .* rep_matrix(baseProbVec, K-1);
  for(n in 1:N) fitted_ref[n] = 1-sum(fitted_nonref[n]);
  fitted = append_col(fitted_ref, fitted_nonref);
  
}

/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 M         = number of covariates
 bg        = established risk (or protective) factors
 tau       = scale parameter

 NM = non-diabetic men
 NW = non-diabetic women
 DM = diabetic men
 DW = diabetic women
 (or any other four groups induced by two binary variables)
*/
// Tomi Peltola, tomi.peltola@aalto.fi

functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] <- sqrt(x[m]);
    }

    return res;
  }

  matrix joint_prior_lp(
        matrix beta_raw, // raw beta parameters
        vector csprime,  // cp and sp
        vector cs_params,
        vector r,        // scales of beta
        matrix V         // eigenvectors of the correlation matrix
        ) {
    matrix[dims(r)[1], 4] beta;

    csprime[1] ~ beta(cs_params[1], cs_params[2]);
    csprime[2] ~ beta(cs_params[3], cs_params[4]);
    col(beta_raw, 1) ~ normal(0.0, 1.0);
    col(beta_raw, 2) ~ normal(0.0, 1.0);
    col(beta_raw, 3) ~ normal(0.0, 1.0);
    col(beta_raw, 4) ~ normal(0.0, 1.0);

    {
      real c;
      real s;
      real lambda;
      vector[4] eigvals;

      c <- 1.0 - 1.0/(1.0 - csprime[1]);
      s <- 1.0 - 1.0/(1.0 - csprime[2]);
      lambda <- sqrt((2.0 * c - 1.0) * (2.0 * s - 1.0) * (2.0 * (c + s) - 1.0) / ((1.0 - 2.0 * (c + s - c * s)) * (c + s - 1.0)));

      eigvals[1] <- 1.0;
      eigvals[2] <- 1.0/sqrt(1.0 - 2.0 * c);
      eigvals[3] <- 1.0/sqrt(1.0 - 2.0 * s);
      eigvals[4] <- 1.0/sqrt(1.0 - 2.0 * (c + s));

      for (m in 1:dims(r)[1]) {
        beta[m] <- (r[m] * lambda) * (V * (eigvals .* (V' * beta_raw[m]')))';
      }
    }

    return beta;
  }

  vector bg_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
  }
}

data {
  int<lower=0> NobsNM;
  int<lower=0> NobsNW;
  int<lower=0> NobsDM;
  int<lower=0> NobsDW;
  int<lower=0> NcenNM;
  int<lower=0> NcenNW;
  int<lower=0> NcenDM;
  int<lower=0> NcenDW;
  int<lower=0> M_bg;
  vector[NobsNM] yobsNM;
  vector[NcenNM] ycenNM;
  matrix[NobsNM, M_bg] Xobs_bgNM;
  matrix[NcenNM, M_bg] Xcen_bgNM;
  vector[NobsNW] yobsNW;
  vector[NcenNW] ycenNW;
  matrix[NobsNW, M_bg] Xobs_bgNW;
  matrix[NcenNW, M_bg] Xcen_bgNW;
  vector[NobsDM] yobsDM;
  vector[NcenDM] ycenDM;
  matrix[NobsDM, M_bg] Xobs_bgDM;
  matrix[NcenDM, M_bg] Xcen_bgDM;
  vector[NobsDW] yobsDW;
  vector[NcenDW] ycenDW;
  matrix[NobsDW, M_bg] Xobs_bgDW;
  matrix[NcenDW, M_bg] Xcen_bgDW;
}

transformed data {
  real<lower=0> tau_mu;
  vector<lower=0>[1] tau_al;
  matrix[4,4] V; // eigenvectors of the correlation matrix

  tau_mu <- 10.0;
  tau_al[1] <- 10.0;

  V[1,1] <-  0.5; V[2,1] <-  0.5; V[3,1] <-  0.5; V[4,1] <- 0.5;
  V[1,2] <- -0.5; V[2,2] <-  0.5; V[3,2] <- -0.5; V[4,2] <- 0.5;
  V[1,3] <- -0.5; V[2,3] <- -0.5; V[3,3] <-  0.5; V[4,3] <- 0.5;
  V[1,4] <-  0.5; V[2,4] <- -0.5; V[3,4] <- -0.5; V[4,4] <- 0.5;
}

parameters {
  vector<lower=0,upper=1>[2] csprime_bg;
  vector<lower=0,upper=1>[2] csprime_al;
  vector<lower=0>[4] cs_params; // a_c & b_c, and a_s & b_s

  real<lower=0> tau_s_bg_raw;
  vector<lower=0>[M_bg] tau_bg_raw;

  matrix[1, 4] alpha_raw;
  matrix[M_bg, 4] beta_bg_raw;

  vector[4] mu;
}

transformed parameters {
  matrix[M_bg,4] beta_bg;
  matrix[1,4] alpha;

  beta_bg <- joint_prior_lp(beta_bg_raw, csprime_bg, cs_params,
                            bg_prior_lp(tau_s_bg_raw, tau_bg_raw), V);

  alpha <- exp(joint_prior_lp(alpha_raw, csprime_al, cs_params, tau_al, V));
}

model {
  yobsNM ~ weibull(alpha[1,1], exp(-(mu[1] + Xobs_bgNM * col(beta_bg, 1))/alpha[1,1]));
  yobsNW ~ weibull(alpha[1,2], exp(-(mu[2] + Xobs_bgNW * col(beta_bg, 2))/alpha[1,2]));
  yobsDM ~ weibull(alpha[1,3], exp(-(mu[3] + Xobs_bgDM * col(beta_bg, 3))/alpha[1,3]));
  yobsDW ~ weibull(alpha[1,4], exp(-(mu[4] + Xobs_bgDW * col(beta_bg, 4))/alpha[1,4]));
  increment_log_prob(weibull_ccdf_log(ycenNM, alpha[1,1], exp(-(mu[1] + Xcen_bgNM * col(beta_bg, 1))/alpha[1,1])));
  increment_log_prob(weibull_ccdf_log(ycenNW, alpha[1,2], exp(-(mu[2] + Xcen_bgNW * col(beta_bg, 2))/alpha[1,2])));
  increment_log_prob(weibull_ccdf_log(ycenDM, alpha[1,3], exp(-(mu[3] + Xcen_bgDM * col(beta_bg, 3))/alpha[1,3])));
  increment_log_prob(weibull_ccdf_log(ycenDW, alpha[1,4], exp(-(mu[4] + Xcen_bgDW * col(beta_bg, 4))/alpha[1,4])));

//    cs_params ~ gamma(0.5, 0.25);
// changed to half-normal
  cs_params ~ normal(0.0, 10.0);

  mu ~ normal(0.0, tau_mu);
}

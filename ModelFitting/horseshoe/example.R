library(rstan)

# data parameters
NobsNM <- 20
NobsNW <- 20
NobsDM <- 20
NobsDW <- 20
NcenNM <- 20
NcenNW <- 20
NcenDM <- 20
NcenDW <- 20
M_bg <- 2
M_biom <- 10
generate_pooled_data <- T

set.seed(123456)

# generate data
bbg <- list()   # established risk factor coeffs
bbiom <- list() # candidate risk factor coeffs
if (generate_pooled_data){
  bbg[[1]] <- rnorm(M_bg)
  bbiom[[1]] <- rnorm(M_biom)/rgamma(M_biom, 5)
  # set half to null to simulate sparsity
  bbiom[[1]][(ceiling(M_biom/2)+1):M_biom] <- 0
  bbg[[2]] <- bbg[[1]]; bbg[[3]] <- bbg[[1]]; bbg[[4]] <- bbg[[1]]
  bbiom[[2]] <- bbiom[[1]]; bbiom[[3]] <- bbiom[[1]]; bbiom[[4]] <- bbiom[[1]]
  alfa <- (1 + runif(1)) * c(1,1,1,1);
} else {
  for (i in 1:4){
    bbg[[i]] <- rnorm(M_bg)
    bbiom[[i]] <- rnorm(M_biom)
    # set half to null to simulate sparsity
    bbiom[[i]][(ceiling(M_biom/2)+1):M_biom] <- 0
  }
  alfa <- 1.0 + runif(4);
}

randn <- function(n, m) matrix(rnorm(n*m), nrow=n, ncol=m)
rand_times <- function(n) 1.2 * runif(n)

data <- list(M_bg = M_bg, M_biom = M_biom)

i <- 0
for (g in c("NM", "NW", "DM", "DW")){
  i <- i + 1
  for (o in c("obs", "cen")){
    names <- c(sprintf("N%s%s", o, g), sprintf("X%s_bg%s", o, g), sprintf("X%s_biom%s", o, g), sprintf("y%s%s", o, g))
    data[[names[1]]] <- get(names[1]);
    data[[names[2]]] <- randn(get(names[1]), M_bg)
    data[[names[3]]] <- randn(get(names[1]), M_biom)
#    data[[names[4]]] <- rand_times(get(names[1]))
    if (o == "obs") data[[names[4]]] <- rweibull(get(names[1]), alfa[i], exp(-(data[[names[2]]] %*% bbg[[i]] + data[[names[3]]] %*% bbiom[[i]])/alfa[i]))
    if (o == "cen") data[[names[4]]] <- runif(get(names[1])) * rweibull(get(names[1]), alfa[i], exp(-(data[[names[2]]] %*% bbg[[i]] + data[[names[3]]] %*% bbiom[[i]])/alfa[i]))
  }
}
data1 <- list(M_bg = M_bg, M_biom = M_biom)
data1$Nobs <- data$NobsNM; data1$Ncen <- data$NcenNM;
data1$yobs <- data$yobsNM; data1$ycen <- data$ycenNM;
data1$Xobs_bg <- data$Xobs_bgNM; data1$Xcen_bg <- data$Xcen_bgNM;
data1$Xobs_biom <- data$Xobs_biomNM; data1$Xcen_biom <- data$Xcen_biomNM;

# initial parameters
nchains <- 2
niter <- 1000
nwarmup <- 500

init <- list()
for (i in 1:nchains){
  init[[i]] <- list(
    csprime_biom = runif(2),
    csprime_bg = runif(2),
    csprime_al = runif(2),
    cs_params = abs(rnorm(4)),
    tau_s_bg_raw = 0.1*abs(rnorm(1)),
    tau_bg_raw = abs(rnorm(M_bg)),
    tau_s1_biom_raw = 0.1*abs(rnorm(1)),
    tau_s2_biom_raw = 0.1*abs(rnorm(1)),
    tau_biom_raw = abs(rnorm(M_biom)),
    tau1_biom_raw = abs(rnorm(M_biom)),
    tau2_biom_raw = abs(rnorm(M_biom)),
    alpha_raw = 0.01*randn(1, 4),
    beta_bg_raw = randn(M_bg, 4),
    beta_biom_raw = randn(M_biom, 4),
    mu = rnorm(4)
  )
}
init1 <- list()
for (i in 1:nchains){
  init1[[i]] <- list(
    tau_s_bg_raw = 0.1*abs(rnorm(1)),
    tau_bg_raw = abs(rnorm(M_bg)),
    tau_s1_biom_raw = 0.1*abs(rnorm(1)),
    tau_s2_biom_raw = 0.1*abs(rnorm(1)),
    tau_biom_raw = abs(rnorm(M_biom)),
    tau1_biom_raw = abs(rnorm(M_biom)),
    tau2_biom_raw = abs(rnorm(M_biom)),
    alpha_raw = 0.01*rnorm(1),
    beta_bg_raw = rnorm(M_bg),
    beta_biom_raw = rnorm(M_biom),
    mu = rnorm(1)
  )
}

# run joint models
fit_bg <- stan(file = "wei_bg_joint.stan", data = data, init = init, chains = nchains, iter = niter, warmup = nwarmup)
fit_hs <- stan(file = "wei_hs_joint.stan", data = data, init = init, chains = nchains, iter = niter, warmup = nwarmup)
fit_lap <- stan(file = "wei_lap_joint.stan", data = data, init = init, chains = nchains, iter = niter, warmup = nwarmup)
fit_gau <- stan(file = "wei_gau_joint.stan", data = data, init = init, chains = nchains, iter = niter, warmup = nwarmup)
# run separate model for NM data
fit_bg1 <- stan(file = "wei_bg.stan", data = data1, init = init1, chains = nchains, iter = niter, warmup = nwarmup)
fit_hs1 <- stan(file = "wei_hs.stan", data = data1, init = init1, chains = nchains, iter = niter, warmup = nwarmup)
fit_lap1 <- stan(file = "wei_lap.stan", data = data1, init = init1, chains = nchains, iter = niter, warmup = nwarmup)
fit_gau1 <- stan(file = "wei_gau.stan", data = data1, init = init1, chains = nchains, iter = niter, warmup = nwarmup)

# rudimentary plot
pars_hs <- extract(fit_hs)
pars_lap <- extract(fit_lap)
pars_gau <- extract(fit_gau)
pars_bg <- extract(fit_bg)
pars_hs1 <- extract(fit_hs1)
pars_lap1 <- extract(fit_lap1)
pars_gau1 <- extract(fit_gau1)
pars_bg1 <- extract(fit_bg1)

plt_bg <- function (x, pars, i, j, col) {
 lines(c(x, x), quantile(pars$beta_bg[,i,j], c(0.05, 0.95)), col=col)
 points(x, median(pars$beta_bg[,i,j]), col=col)
}
plt_biom <- function (x, pars, i, j, col) {
 lines(c(x, x), quantile(pars$beta_biom[,i,j], c(0.05, 0.95)), col=col)
 points(x, median(pars$beta_biom[,i,j]), col=col)
}
plt_bg1 <- function (x, pars, i, col) {
 lines(c(x, x), quantile(pars$beta_bg[,i], c(0.05, 0.95)), col=col)
 points(x, median(pars$beta_bg[,i]), col=col)
}
plt_biom1 <- function (x, pars, i, col) {
 lines(c(x, x), quantile(pars$beta_biom[,i], c(0.05, 0.95)), col=col)
 points(x, median(pars$beta_biom[,i]), col=col)
}


ofs <- c(-0.3, -0.1, 0.1, 0.3)
plot(c(0, M_bg+M_biom+1), c(0, 0), xlim=c(0, M_bg+M_biom+1), ylim=c(-2.5,2.5), type='l')
for (i in 1:M_bg){
  for (j in 1:4){
    plt_bg(i+ofs[j]-0.05, pars_hs, i, j, 'blue')
    plt_bg(i+ofs[j], pars_lap, i, j, 'black')
    plt_bg(i+ofs[j]+0.05, pars_gau, i, j, 'red')
    plt_bg(i+ofs[j]+0.1, pars_bg, i, j, 'magenta')
  }
  j <- 1
  plt_bg1(i+ofs[j]-0.03, pars_hs1, i, 'lightblue')
  plt_bg1(i+ofs[j]+0.02, pars_lap1, i, 'lightgray')
  plt_bg1(i+ofs[j]+0.07, pars_gau1, i, 'lightpink1')
  plt_bg1(i+ofs[j]+0.12, pars_bg1, i, 'lightsalmon')
}
for (i in 1:M_biom){
  for (j in 1:4){
    plt_biom(i+ofs[j]+M_bg-0.05, pars_hs, i, j, 'blue')
    plt_biom(i+ofs[j]+M_bg, pars_lap, i, j, 'black')
    plt_biom(i+ofs[j]+M_bg+0.05, pars_gau, i, j, 'red')
  }
  j <- 1
  plt_biom1(i+ofs[j]+M_bg-0.03, pars_hs1, i, 'lightblue')
  plt_biom1(i+ofs[j]+M_bg+0.02, pars_lap1, i, 'lightgray')
  plt_biom1(i+ofs[j]+M_bg+0.07, pars_gau1, i, 'lightpink1')
}

